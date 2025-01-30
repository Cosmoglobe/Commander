	Program FTB_HISTO

C------------------------------------------------------------------------
C
C    PURPOSE: Plots a histogram of raw interferogram counts.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: R. Isaacman
C            ARC
C            September 1, 1988
C
C    INVOCATION: HISTO
C
C    INPUT PARAMETERS: None
C
C    OUTPUT PARAMETERS:
C
C    SUBROUTINES CALLED:
C	FTB_accum_HISTO
C
C    COMMON VARIABLES USED:
C
C    INCLUDE FILES:
C	$SSDEF
C	CT Lib
C
C----------------------------------------------------------------------
C    CHANGE LOG:
C
C	Shirley M. Read, STX, September 7, 1988.
C		Converted program to use the PLT package instead of
C	        Vecplt in response to a decision of the FIRAS Team.
C               SPRs 1566 and 1763 requested that all FIRAS facilities
C               should tell of successful completion and signal errors
C	        via $status for batch processing. Also added interface
C		with the FUT_Error, condition handler.
C	RAShafer, Code 685, October 5 1988.
C		Change read form to FTB_ACCUM_HISTO to allow more than
C		100 histograms to be used.
C
C	Version 4.2.1 11/08/88, SPR 2737, Shirley M. Read, STX
C		Replace call to XQUEST with standard FIRAS user interface.
C
C
C       Version 4.4   August 7, 1989, Jon Smid, STX
C               SPR 3663, Implementation of the noninteractive mode
C
C       Version 4.5   1989 Sep 25,  SPR 4517,  Fred Shuman, STX.
C               Clear yvals() array each time a new timerange is requested.
C               FTB_HISTO and FTB_ACCUM_HISTO.
C       SPR 4900, November 3, 1989, Aliza R. Panitz, STX
C               Get start/stop times from FTB_ACCUM_HISTO and put on 
C               plot labels.
c
c       SPR 4797, Nov 8, 1989, D. Bouler, STX
c               Signal actual value in istat at end of program rather than
c               just ftb_aberr.
c
c	SER 5726, Mar 6, 1990, S. Alexander, STX
c		Add capability to pass in a PLT command file.
c
C----------------------------------------------------------------------

	implicit none

C	FTB messages

	external	ftb_normal
	external	ftb_aberr
	external	ftb_ctinit
	external        fut_normal
	external        cli$_present
	external        cli$_defaulted

C	Condition Handler

	external        fut_error

C	Include Files

	include '($ssdef)'
	Include 'CT$Library:CTUser.Inc'


C	Local Declarations

	real		*4	xvals(4096)
	real		*4	yvals(4096)


	integer         *4      cli$present
	integer         *4      cli$get_value
	integer         *4      status
	integer		*4	ichan
	integer		*4	istat
	integer		*4	n
	integer		*4	j
	integer		*4	k
	integer		*4	nbin
	integer		*4	istart
	integer		*4	iend
	integer		*4	nraw
	integer		*4	numrec
	logical		*1	gotchan / .false. /
	logical		*4	more_time
	character	*4	ans
        character       *14     jstart
        character       *14     jstop
	Integer * 2	ct_stat(20)	!CT return status
	Integer         *4      len

C	Information for PLT call.

	real		*4      y(4096,3)     ! Vectors to be plotted
	real		*4	ymin          ! Minimum Y value for plot
	real		*4	ymax          ! Maximum Y value for plot
	real		*4	delta	      ! Difference of ymin and ymax
	real		*4	ypltmin	      ! Min used in the plot
	real		*4	ypltmax	      ! Max used in the plot
	integer		*4	mx1 /4096/    ! Actual 1st dimension of Y array
	integer		*4      npts          ! Number of points to plot
	integer		*4	nveci /2/     ! Number of vectors to be plotted
	integer		*4	iery(2) /0,1/ ! Flags for plot errors: 0=no err
	integer		*4	ier	      ! Error flag
	integer		*4	ncmd          ! Number of commands
	character	*150    cmd(50)       ! Plot commands
	data	        cmd(2)  / 'LA X IFG Value (Raw Counts): ' /
	data	        cmd(3)  / 'LA Y Number of Occurrences' /
	character	*28     title / 'LA T Histogram of Raw IFGs: '/
	character	*29     xlabel/ 'LA X IFG Value (Raw Counts): ' /
	character	*11     lchan / 'Channel = ' /
	character	*2      chan_ch(4) / 'RH', 'RL', 'LH', 'LL' /
	character       *2      chan
	character       *50     plotfile
	character	*64	plt_com_file
	character       *3      input

	logical         *1      inter


C	Functions

	integer		*4	ftb_accum_histo
	integer		*4      fut_erase
	integer         *4      str$upcase

C
C     Establish condition handler.
C
	call lib$establish ( fut_error )
c
c Get the channel
c
C
C Get inputs.
C
	If ( cli$present('interactive').eq.%loc(cli$_present) .or.
     &       cli$present('interactive').eq.%loc(cli$_defaulted) ) then
	   Inter = .true.
	else
	   Inter = .false.
	endif

	If ( .not. inter ) then
	   status = CLI$Get_Value('Input',Input,Len)
	   If  ( input .eq. 'RAW' )  then
	      title = 'LA T Histogram of Raw IFGs: '
	   else if (input .eq. 'FPP') then
	      title = 'LA T Histogram of Raw FPP IFGs:'
	   else
	      title = 'LA T Histogram of Raw FDQ IFGs:'
	   endif
	   status = CLI$Get_Value('Channel',Chan,Len)
	   If ( chan .eq.'RH' ) then
	      ichan = 1
	   else if ( chan .eq.'RL' ) then
	      ichan = 2
	   else if ( chan .eq.'LH' ) then
	      ichan = 3
	   else
	      ichan = 4
	   endif
	endif

c
c  Interactive or noninteractive branch
c
	If ( inter ) then
	   do while ( .not. gotchan )
	      type 23
23            format (' Enter channel number: '$)
	      accept *, ichan

	      if ((ichan .ge. 1) .and. (ichan .le. 4)) then
!	data	        cmd(2)  / 'LA X IFG Value (Raw Counts) ' /
!	         cmd(1) = title // lchan // chan_ch(ichan)
                 cmd(2) = xlabel // lchan // chan_ch(ichan)
	         gotchan = .true.
	      else
	         type 24
24               format (' Valid channels are 1 to 4')
	      endif
	   enddo

	else
!!	   cmd (1) = title // lchan // chan
           cmd(2) = xlabel // lchan // chan
	endif
c
c   Initialize COBETRIEVE
c
	call ct_init(ct_stat)
	if (ct_stat(1) .ne. ctp_normal) then
	   call lib$signal(ftb_ctinit, %val(1), %val(ct_stat(1)))
	   call lib$signal(ftb_aberr)
	   call exit(ss$_abort)
	end if

	more_time = .true.
      	do while(more_time)

	   do j=1,4096
	      yvals(j) = 0.
	      xvals(j) = j - 2049
	   enddo

	   istat = ftb_accum_histo (ichan, numrec, yvals, inter,
     ^                              jstart,jstop)

  	   if (istat .eq. %loc(ftb_normal)) then
C
C Find nonzero part of histogram
C
	      istart = 1
	      do while (yvals(istart) .lt. 0.1)
	         istart = istart + 1
	      enddo

C  Convert istart to actual sample No.

	      istart = istart - 10 - 2049

	      iend = 4096
	      do while (yvals(iend) .lt. 0.1)
	         iend = iend - 1
	      enddo
	      iend = iend + 10 - 2049

	      npts = 4096
	      nraw = 0

C Calculate total number of samples

	      do j=1,npts
	         nraw = nraw + yvals(j)	
	      enddo

C Load the PLT y matrix.

	      do j = 1,npts
	         y(j,1) = xvals(j)
	         y(j,2) = yvals(j)
	         if(yvals(j) .eq. 0) then
	            y(j,3) = .5
	         else
	            y(j,3) = sqrt(yvals(j))
	         end if
	      enddo

C Determine y min and y max.

	      ymin = y(1,2)
	      ymax = y(1,2)
	      do j = 2, npts
	         if ( y(j,2) .lt. ymin ) ymin = y(j,2)
	         if ( y(j,2) .gt. ymax ) ymax = y(j,2)
	      enddo
	      delta = 1.05 * (ymax-ymin)
	      if (delta .gt. 0) then
	         ypltmax = ymin + delta
	         ypltmin = ymax - delta
	      else
	         ypltmax = ymax + 0.1
	         ypltmin = ymin - 0.1
	      endif
	      write (cmd(4),'(a, 2(1x,e14.4))') 'R Y ', ypltmin, ypltmax
       	      write (cmd(5),'(a, 2i6)') 'R X ',istart,iend
	      ncmd = 5

	      If  ( .not. inter ) then
	         status =  Cli$get_value ('PLOTDEVICE', plotfile,len )
	         ncmd = ncmd + 1
	         cmd(ncmd) = 'D '// plotfile
	      endif

	      if (cli$present('PLTFILE') .eq. %loc(cli$_present)) then
	         status = cli$get_value('PLTFILE',plt_com_file,len)
	         ncmd = ncmd + 1
	         cmd(ncmd) = '@' // plt_com_file
	      end if

	      If  ( inter ) then
	         type *, 'Number of records =',numrec,'  Number of points =',nraw
	      else
	         ncmd = ncmd + 2
	         cmd(ncmd-1) = 'P'
	         cmd(ncmd) = 'Q'
	      endif

              cmd(1) = title // jstart // ' To ' // jstop
C Invoke PLT
C!!! scratch memo area
C! the titles and labels are in the CMD array...
C!!! end scratch

	      call plt(y,iery,mx1,npts,nveci,cmd,ncmd,ier)

	      istat = fut_erase()
	      if (istat .eq. %loc(fut_normal)) istat = %loc(ftb_normal)
	   end if

C Inquire if user wants another time range. Replace Xquest call with
C standard FORTRAN interface calls.

	   If ( inter ) then
	      type 25
25            format (' Another time range? ([Y]/N) '$)
	      accept 10, ans
10            format (a)
	      status = str$upcase (ans,ans)
	      if(ans(1:1) .eq. 'N') more_time = .false.
	   else
	      more_time = .false.
	   endif

	end do


C	Exit the program.

	if (istat .eq. %loc(ftb_normal)) then
	   call lib$signal(ftb_normal)
	   call exit(ss$_normal)
	else
           call lib$signal (%val(istat))
	   call lib$signal(ftb_aberr)
	   call exit(ss$_abort)
	endif

	end
