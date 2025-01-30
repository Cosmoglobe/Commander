	program ftb_mtm_jitter

c-------------------------------------------------------------------------
c
c  Program FTB_MTM_JITTER
c
c  Analysis program for display of MTM jitter data taken in
c  microprocessor mode D1.
c
c  Written 12 August 1986
c  Rich Isaacman
c  Applied Research Corp.
c
c  Include files used:
c    $ssdef
c
c  Modified:
c
c  Reid Wilson,  STX Inc.,  21 November 1987
c
c  SPR 1346, Fetch data from the archive.  R. Kummerer, December 17, 1987.
c
c  SPR 1351, Convert to new IMSL. R. Kummerer, April 14, 1988.
c
c
c  Modified 15 Aug 1988 by RBI to accept mode 5 (Seismo) data
c
c  Modified August 25, 1988
c  Shirley M. Read
c  STX Inc.
c    Reason: SPRs 1566 and 1763 requested that all FIRAS facilities
c	     should tell of successful completion and signal errors
c	     via $status for batch processing. Also added the interface
c	     to FUT_Error, the condition handler.
c            Configured to meet CSDR standards.
c
c  Modified September 26, 1988
c       Shirley M. Read
c       STX Inc.
c       Reason: Added capability to produce automatic hardcopies, and include
c	        additional information in plot labels.(SPR 2430)
c
c  Modified October 5, 1988
c       Shirley M. Read, STX
c       Reason: Corrected the inclusion of the FTB_Labels text file in the
c               subroutine. (SPR 2558)
c
c  Modified May 2, 1989
c       R. Kummerer, STX
c       Reason: Use FUT_VECPLT instead of FTB_VECPLT. (SPR 3769)
c
c  Modified March 6, 1990
c	S. Alexander, STX
c	Reason: Add capability to pass in a PLT command file. (SER 5726)
c
c  F. Shuman, STX, 1990 Mar 09:  SPRs 5396, 5410 -- recent Jitter runs
c     disagree with old ones.  Module FTB_PlotSpec.for
c
c  Modified 15 June 1992, L. Rosen, Hughes STX. SER 2865, SPR 9765.
c     Add plot of FFT of average of samples.
c-------------------------------------------------------------------------

        implicit none

c	FTB messages

	external		fut_normal
	external		ftb_normal
	external		ftb_aberr
	external		ftb_ctinit
	external		ftb_ctopen
	external		ftb_ctread
	external		ftb_ctclos

c	Condition Handler

	external        fut_error

c	Include files.

	include		'ct$library:ctuser.inc'
	include		'(fut_vecplt_labels)'
	include		'(fut_params)'
	include		'(upm_stat_msg)'
	include		'($ssdef)'

c	Local Declarations

	integer		*2      ct_nr
	integer		*2      ctstat(20)
	integer		*4	ibin
	integer		*4	ios
	integer		*4	status

	integer		*4      start_time(2)
	integer		*4      end_time(2)
	character	*14     starting_time
	character	*14	ending_time
	logical		*1      prompt_times
	logical		*1      seeall
	integer		*4      lstart
	integer		*4      lstop
	character	*14     jstart_default
	character	*14     jstop_default

	integer		*4	plt_com
	character	*64	plt_com_file

        real		*4	grand_spec(512)
        real		*4	histo(256)
        real		*4	time(1024)
        real		*4	x(1024)
        real		*4	xbuff(1024)
        real		*4	timebuff(1024)
	real		*4	samples(1024)
	complex		*8	zspec(1024)		! complex spectra
	complex		*8	avg_zspec(512)		! sum and average z spec
	real		*4	apod			! apodization function
	real		*4	fft_avg(512)		! new fft of avg.
	real		*4	average_time		! new
        real		*4	sumtime
        real		*4	sumsq
        real		*4	avg_time
        real		*4	rms_time
        real		*4	grand_avg
        real		*4	grand_rms
	character	*100	text
	character	*100	filename
	character	*16	time_tag
        integer		*4	nrec
        integer		*4	ittag
        integer		*4	j
        integer		*4	len
        integer		*4	kstart
	integer		*4	chan
        integer		*4	isample(1024)
        integer		*4	n
        integer		*4	m
        integer		*4	npts
        integer		*4	k
	integer		*4	mode_d
	integer		*2	map
	byte			bmap(2)
	byte			yn
	logical		*1	normal/.true./
	integer         *4      plotmode       !Screen display or auto hardcopy

	common /worksp/ rwksp	! Work space initialization for FFTCF
	real		*4	rwksp(6179)

c	Functions

	integer		*4	FUT_Vecplt
	integer		*4      fut_timerange
	integer		*4      upm_get_value
	integer		*4      upm_present
	integer		*4	ftb_mtm_sumplot

	dictionary 'nfs_emf'
	record /nfs_emf/ emf_rec

	data jstart_default/'85001000000000'/
	data jstop_default/'99365235959990'/

	equivalence (map, bmap(1))


c
c Establish condition handler.
c
	call Lib$Establish( FUT_Error )

c
c Initialization.

	call iwkin(6179)	! Work space initialization for FFTCF

	title = 'MTM Jitter: YY:DDD:HH:MM:SS '
	oxlabl = ' '
	oylabl = ' '

c	call ustart

c
c Get the JSTART and JSTOP times if running from an automation command
c file.
c
	starting_time = jstart_default
	ending_time = jstop_default
	prompt_times = .true.

	if (upm_present ( 'JSTART' ) .eq. UPM_Pres) then
           ios = upm_get_value ( 'JSTART', starting_time, len )
           if(ios .eq. ss$_normal)then
	      lstart = index(starting_time,' ') - 1
	      if(lstart .eq. -1)lstart=14
	      starting_time = starting_time(1:lstart) //
     .				jstart_default(lstart+1:)
           end if
	   prompt_times = .false.
	end if

	if (upm_present ( 'JSTOP' ) .eq. UPM_Pres) then
           ios = upm_get_value ( 'JSTOP', ending_time, len )
           if(ios .eq. ss$_normal)then
	      lstop = index(ending_time,' ') - 1
	      if(lstop .eq. -1)lstop=14
	      ending_time = ending_time(1:lstop) //
     .				jstop_default(lstop+1:)
           end if
	   prompt_times = .false.
	end if

	time_range = starting_time//';'//ending_time//';'

        if (prompt_times) then
	   ios = fut_timerange ( starting_time,
     .			start_time, ending_time, end_time )
	   time_range = starting_time//';'//ending_time//';'
	end if
c
c Retrieve PLT command file if specified.
c
	plt_com = fac_not_present

	if (upm_present ( 'PLTFILE' ) .eq. upm_pres) then
	   plt_com = fac_present
	   status = upm_get_value( 'PLTFILE', plt_com_file, len)
	end if
c
c Fetch a channel number.
c
	chan = 0

	do while (chan .lt. 1 .or. chan .gt. 4)
	   type 30
30	   format(/,' Enter a channel number: ', $)
	   accept 45, chan
45	   format (i)
	end do
c
c Prompt for full or abbreviated display.
c
	type 60
60	format (/' Do you wish to see every individual sweep (Y/[N])? '$)
	accept 65, yn
65	format (a1)
	seeall = .false.
	if (yn.eq.'y' .or. yn.eq.'Y') seeall = .true.
c
c  Prompt for "regular" jitter or seismographic mode data
c
	type 72
72	format (/' Do you want [J]itter (default) or [S]eismo data? '$)
	accept 65, yn
	mode_d = 1
	if (yn.eq.'s' .or. yn.eq.'S') mode_d = 5
	if (mode_d.eq.1)then
	write(title2,900)chan
900	format('Channel: ',i1,' Mode: Jitter')
	else
	write(title2,901)chan
901	format('Channel: ',i1,' Mode: Seismo')
	endif
c
c  Prompt for screen display or automatic hardcopy.
c
	type 80
80	format (/' Do you want [C]rt display (default) or ',
	1 '[A]utomatic hardcopy? '$)
	accept 65, yn
	plotmode = 1
	if (yn.eq.'a' .or. yn.eq.'A') plotmode = 2
c
c Initialize Cobetrieve.
c
	call ct_init(ctstat)

	if (ctstat(1) .eq. ctp_normal) then

	   if (chan .eq. 1) then
	      call ct_open_arcv(,ct_nr,ctu_$firas,ctu_$fir_ed1,ctu_$r,
     .			      ctstat,%ref(time_range),)
	   else if (chan .eq. 2) then
	      call ct_open_arcv(,ct_nr,ctu_$firas,ctu_$fir_ed2,ctu_$r,
     .			      ctstat,%ref(time_range),)
	   else if (chan .eq. 3) then
	      call ct_open_arcv(,ct_nr,ctu_$firas,ctu_$fir_ed3,ctu_$r,
     .			      ctstat,%ref(time_range),)
	   else if (chan .eq. 4) then
	      call ct_open_arcv(,ct_nr,ctu_$firas,ctu_$fir_ed4,ctu_$r,
     .			      ctstat,%ref(time_range),)
	   end if

	   if (ctstat(1) .eq. ctp_normal) then

c
c Look for diagnostic mode data
c
	      nrec = 0
	      do j=1,256
	         histo(j) = 0.
	      enddo

	      call ct_read_arcv ( , ct_nr, emf_rec, ctstat )

	      do while ((ctstat(1) .eq. ctp_normal) .and. (normal))

		 if (emf_rec.sci_head.sc_head1a .eq. mode_d .and.
     .		     emf_rec.sci_head.sc_head1b .eq. 'D') then

	            title(13:28) = emf_rec.ct_head.gmt

	            nrec = nrec + 1

c
c		Decode the timing information.
c
	            do j=1,512
	               map = emf_rec.sc_eng_data(j)
	               k = j*2
	               isample(k-1) = bmap(2)
	               isample(k) = bmap(1)
	            enddo

	            do j=1,1024
	               if (isample(j) .lt. 0) then
	                  isample(j) = isample(j) + 256
	               end if
	            enddo

c
c		Convert times to microsec, then calculate means and RMS.
c		Accumulate times in histogram buffer "histo". Save times
c		in 2 places since "time" arry gets overwritten by plotspec.
c
	            sumtime = 0.
	            sumsq = 0.

	            do n=1,1024
	               time(n) = 10.*isample(n) + 45.
		       samples(n) = time(n)
		       ibin = isample(n) + 1
		       histo(ibin) = histo(ibin) + 1
	               x(n) = n
	               sumtime = sumtime + time(n)
	               sumsq = sumsq + time(n)*time(n)
	            enddo

	            avg_time = sumtime/1024.
		    average_time = average_time + avg_time
	            rms_time =  amax1((sumsq/1024. - avg_time*avg_time),0.)
		    rms_time = sqrt (rms_time)
		    grand_avg = (grand_avg*(nrec-1) + avg_time)/nrec
		    grand_rms = sqrt ((grand_rms**2*(nrec-1) +
     1					rms_time**2)/nrec)
c
c		Accumulate average of FFT (spectrum)
c
		    call ftb_plotspec (time, avg_time, timebuff,
     1				       xbuff, .false., plotmode,
     2				       plt_com, plt_com_file)
		    do j=1,512
		       grand_spec(j) = (grand_spec(j)**2*(nrec-1) +
     1					time(j)*time(j))/nrec
		       grand_spec(j) = sqrt(grand_spec(j))
		    enddo
		    do j=1,1024
		       apod = 2./3.*( 1. + cos((j-512.5)*3.1415927/511.5) )
		       zspec(j) = cmplx (samples(j)-avg_time,0.) * apod
		    enddo
		    call fftcf (1024, zspec, zspec)
		    do j=1,512
		       avg_zspec(j) = avg_zspec(j) + zspec(j)
		    enddo
c
c		Display results for this sweep only if
c		"seeall" flag has been set. Then loop back
c		for next record.
c
		    if (seeall) then

	               do m=1,1024
	                  xbuff(m) = x(m)
	                  timebuff(m) = samples(m)
	               enddo
		       xlabl = 'Sample Number'
		       ylabl = 'Sample Interval (\gmsec)'
		       write(oxlabl,'(a,f6.1,a,f6.1,a)')
     &			'Mean Interval:',avg_time,' +/- ',rms_time,
     &			' (rms) \gmsec'

	               npts = 1024

	               status = fut_vecplt (xbuff,timebuff,npts,
     .			 plotmode,0,0,0,0,0,0,0,0, plt_com, plt_com_file)
		       oxlabl = ' '
		       if ( status .ne. %loc(FUT_Normal)) then
			 normal = .false.
		       else
	                 call ftb_plotspec (samples, avg_time, timebuff,
     1					    xbuff, .true., plotmode,
     2					    plt_com, plt_com_file)

	                 type 100, nrec, emf_rec.ct_head.gmt, avg_time,
     .				rms_time
100	                 format (24(/),27x,'Record number:',i4,//20x,'Time: ',
     .			   a14//,20x,'Average sample:',f6.1,
     .	                   ' microsec'//20x,'RMS jitter:',f6.1,
     .			   ' microsec'/////////)

		       endif
		    endif

		 end if
	         if (normal) call ct_read_arcv(,ct_nr,emf_rec,ctstat)

	      end do
	      if (normal) then
	        if (ctstat(1) .ne. ctp_endoffile) then
	          call lib$signal(ftb_ctread,%val(1),%val(ctstat(1)))
		  normal = .false.
	        else if (nrec .eq. 0) then
		  type 200
200		  format (/, ' No data found for MTM Jitter.', /)
	        else
		   average_time = average_time / nrec
		   do j = 1,512
		      avg_zspec(j) = avg_zspec(j) / nrec
		      fft_avg(j) = sqrt (real (avg_zspec(j) * 
     1		         conjg (avg_zspec(j)))) / 512.
		   enddo
c
c  Make plots
c
		   status = ftb_mtm_sumplot (histo, grand_spec, grand_avg,
     1		      grand_rms, nrec, timebuff, xbuff, plotmode, plt_com,
     2		      plt_com_file, fft_avg, average_time)
	        end if
	      endif

	   else
	      call lib$signal(ftb_ctopen,%val(1),%val(ctstat(1)))
	      normal = .false.
	   end if

	else
	   call lib$signal(ftb_ctinit,%val(1),%val(ctstat(1)))
	   normal = .false.
	end if

	call ct_close_arcv(,ct_nr,ctstat)

	if (ctstat(1) .ne. ctp_normal) then
	   call lib$signal(ftb_ctclos,%val(1),%val(ctstat(1)))
	   normal = .false.
	end if

c       Exit the program.

	if (normal) then
	  call lib$signal(ftb_normal)
	  call exit(ss$_normal)
	else
	  call lib$signal(ftb_aberr)
	  call exit(ss$_abort)
	end if

	end
