	subroutine fsd_glitchmap_data(sci_rec,s_flag,number,size,lun,cindex,
	1                             config,new_segment,stat)

C---------------------------------------------------------------------   
C
C       Version 4.2.1 12/02/88, SPR 2310, Quoc Chung, STX
C		Bring FSD up to standard ,error message trapping,
C       	and exit status.
C	Version 4.4.1 09/15/89, SER 4567, R. Kummerer, STX
C		Use PLT instead of TEMPLATE for graphics.
C	Version 4.4.1 10/06/89, SPR 4137, R. Kummerer, STX
C		Remove /OPTION as it does not and never will do
C		anything.
C
C       Version 7.3.1 12/10/90, SPR 5866, N. Gonzales, STX
C               Change IFG sampling rate from 684.0 to 681.43
C
C	SPR 9847, Update FSD to use new FEX_Samprate reference file using
C               cct_get_config.  N. Gonzales, Hughes STX / 1992 August 5.
C---------------------------------------------------------------------   

        implicit	none

	include		'(fsd_glitchmap)'
	include         '(cct_get_config)'
	include         '($ssdef)'
c
c Functions
c    
	integer   *  4      cct_get_config_tod
c
c Records
c
	Dictionary 'fex_samprate'
	Structure /config_data/
	   Record /fex_samprate/fex_samprate
	Endstructure

	Record /config_data/ config
	Record /config_status/ stat(1)

	Integer   * 4       number   	      ! number of data sets
	Integer   * 4       size              ! size of data sets in bytes
	Integer   * 4       lun(1)	      ! logical unit numbers
	Integer   * 4       cindex(1)         ! initial cache pointers
	Logical   * 1       new_segment(1)    ! flag for new segments
	
	character * 14      start_time
	character * 14      stop_time
	character *  2      scan_mode
	integer   *  4      status
	integer   *  4      speed
	integer   *  4      length
	integer   *  4      ngr,nsw,ibit
	integer   *  4      iglitch,nglitch,igmap
	integer   *  4      ishift
	integer   *  4      fake
	integer   *  4      uP
	integer   *  4      gain_code
	integer   *  4      xpos
	integer   *  2      gltch_w(32)
	integer   *  2      gg(64)
	integer   *  4      jtime(2)

	real      *  4      rifg(512)
	real      *  4      hist(512)
	real      *  4      samp_rate(2)
        real      *  4      fact

	character * 14      gmt_time

	byte                gltch(64)

	logical   *  4      i_flag
	logical   *  4      s_flag
	logical   *  4      glitch(512)
	logical   *  4      scan_flag

        external  fsd_maxifg
        external  fsd_endtime
	external  fsd_getconfigerr

	equivalence(gltch,gltch_w)
c
c Extract instrument mode data from the raw science record.
c
	ngr = sci_rec.sci_head.sc_head9
 	nsw = sci_rec.sci_head.sc_head11
	speed = sci_rec.sci_head.mtm_speed
	length = sci_rec.sci_head.mtm_length
	fake = sci_rec.dq_data.fake
	xpos = sci_rec.dq_data.xcal_pos
	gain_code = sci_rec.sci_head.gain
	uP = sci_rec.sci_head.sc_head1a
	start_time = sci_rec.ct_head.gmt

	Call ct_gmt_to_binary (start_time, jtime)
c
c Access FEX_Samprate reference dataset using CCT_Get_Config. 
c
	status = cct_get_config_tod (jtime, number, size, lun, cindex,
	1                          config, new_segment, stat)

	if (.not. status) then
	    call lib$signal(fsd_getconfigerr,%val(1),%val(status))
	endif

	samp_rate(1) = config.fex_samprate.sampling_rate
	samp_rate(2) = 512.

	if (speed .eq. 0) then
	   scan_mode(2:2)='S'
	else
	   scan_mode(2:2)='F'
	end if

	if (length .eq. 0) then
           scan_mode(1:1)='S'
	else
	   scan_mode(1:1)='L'
	end if
	
	scan_flag=.false.
	do i=1,n_scan
	   if (nscan(i) .eq. scan_mode) scan_flag=.true.	
	enddo
	if (.not. scan_flag) return

	if (num+1 .eq. 1) then
	   scan = scan_mode
	   mtm_speed = speed
	   mtm_length = length
	   ngroup = ngr
	   nsweep = nsw
	   fakeit = fake
	   xcal_pos = xpos
	   upmode = uP
	   gain = fac_gains(gain_code)
	   fact = samp_rate(fakeit+1)/ngroup/nsweep/512
	endif

	if ((scan .ne. scan_mode) .or. (fake .ne. fakeit)) then
	   s_flag=.true.
	   more=.false.
	   print *,'New scan mode.'
	else
	   s_flag=.false.
	endif

        if (.not. s_flag) then
	   num=num+1
	   gmt(num)=sci_rec.ct_head.gmt
	   degl_thre(num)=sci_rec.sci_head.sc_head14
	   n_over(num)=sci_rec.sci_head.sc_head22
	   nglitch = sci_rec.sci_head.sc_head21
	   avg_g(num) = nglitch*fact
c
c Obtain individual glitchmaps
c
	   do j=1,32
	      gltch_w(j)=sci_rec.ifg_data.gltch(j)
	   enddo

	   do i=1,64
	      gg(i)=gltch(i)
	   enddo

	   do i=0,511,8
	      do ibit = 0, 7
		 igmap = i/8 + 1
		 iglitch = i + (8-ibit)
		 glitch(iglitch) = gg(igmap)
		 gg(igmap) = ishft(gg(igmap),-1)
	      enddo
	   enddo

	   do i=1,512
	      hist(i)=0
	      if (glitch(i)) then
		 hist(i)=hist(i)+1
    	      endif
	      his(i,num)=hist(i)
	      g_hist(i)=g_hist(i)+hist(i)
	   enddo
c
c Obtain individual IFGs
c
	   if (ifg_flag) then
	      call lib$movc5(0,,0,2048,rifg)
  	      do i=1,512
	         rifg(i)=floati(sci_rec.ifg_data.ifg(i))
	         sci_ifg(i,num)=rifg(i)
	         g_ifg(i)=g_ifg(i)+rifg(i)
	      enddo

	   endif !ifg_flag

	   if (num .ge. 50) then
	      more=.false.
              call lib$signal(fsd_maxifg,%val(1),%val(num))
	      grp_max=.true.
	   else
	      grp_max=.false.
	   endif

	endif

	if (eof) then
	   call lib$signal(fsd_endtime,%val(1),%val(chan_id))
	   more=.false.
	endif

	if ((num .gt. 0) .and. (.not. more)) then
	   if (ifg_flag) then
	      do i=1,512
	         g_ifg(i)=g_ifg(i)/num
	      enddo
	   endif
	endif

	return
	end
