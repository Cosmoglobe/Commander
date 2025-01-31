
      integer*4 function fss_read_skymap(lu_input, sky_buff, sci_recs,
     .                                   sky_recs, closeold, minimum,
     .                                   bins, num_bins, mpmodes, plist, 
     .                                   inst_qual, attit_qual, oldest, 
     .                                   lu_rep, chan, start_chan, stop_chan,
     .                                   chan_int,temptol,sci_times,max_dihed,
     .                                   report, mtmlist,galaxy,minsun,maxsun,
     .                                   minmoon,maxmoon,minearth,maxearth)

C------------------------------------------------------------------------
C    PURPOSE: Read in input skymap. Filter records as per qualifiers.
C             Return array of records and number of records.
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C
C    INVOCATION: status=fss_read_skymap(lu_input, sky_buff, sci_recs,
C     .                                 sky_recs, closeold, minimum,
C     .                                 bins, num_bins, mpmodes, plist, 
C     .                                 inst_qual, attit_qual, oldest, 
C     .                                 lu_rep, chan, start_chan, stop_chan, 
C     .                                 chan_int,temptol,sci_times,max_dihed,
C     .                                 report, mtmlist,galaxy,minsun,maxsun,
C     .                                 minmoon, maxmoon, minearth, maxearth)
C
C    INPUT PARAMETERS:       I*4   lu_input     The input skymap LU
C                            I*4   sci_recs     Science records read
C                            I*4   closeold     Pull in uncoadded groups ?
C                            I*4   minimum      Minimum size for coadd group
C                            R*4   bins(*)      Temp bin midpoints
C                            I*4   num_bins     Number of temp bins
C                            I*4   mpmodes(*)   Microprocessor modes
C                            i*4   plist(*)     Pixels to process
C                            I*4   inst_qual    Instrument quality
C                            I*4   attit_qual   Attitude quality
C                            c*14  oldest       Earliest acceptable time
C                            I*4   lu_rep       Report file LU
C                            I*4   chan         Channel being processed
C                            I*4   start_chan   Start channel
C                            I*4   stop_chan    Stop channel
C                            I*4   chan_int     Channel interval
C                            R*4   temptol      Temperature tolerance
C                            i*4   sci_times(2,4,fss_max_buff) Science times
C                            r*4   max_dihed    Max dihedral temp
C                            i*4   report       Was report requested?
C                            i*4   mtmlist(4)   MTM modes
C
C			     R*4   Galaxy       Exclude recs with lat < Galaxy
C			     R*4   minsun       Minimum sun angle
C			     R*4   maxsun       Maximum sun angle
C			     R*4   minmoon      Minimum moon angle
C			     R*4   maxmoon      Maximum moon angle
C			     R*4   minearth     Minimum earth angle
C			     R*4   maxearth     Maximum earth angle
C
C
C    OUTPUT PARAMETERS:      array sky_buff     Input skymap records
C                            I*4   sky_recs     Number of input sky records
C 
C    SUBROUTINES CALLED:     csa_read_pixels
C
C    COMMON VARIABLES USED:  None.
C 
C    INCLUDE FILES:          include       '(csa_pixel_input_rec)'
C                            include       '(csa_pixel_output_rec)'
C                            include       '(fss_include)'
C
C----------------------------------------------------------------------
c                      PDL for FSS_READ_SKYMAPS
c                      D. Bouler, JAN, 1990
c
c
c  initialize variables
c
c  initialize pixel list to be read
c
c  call csa_read_pixels to read portion of input skymap into input buffer
c
c  check for same pixelization on command line and in input skymap;
c  if not same, halt with message
c
c  check skymap records for consistency; records must be:
c    1.  In the given temperature bins
c    2.  In the microprocessor mode list
c    3.  In the pixel list
c    4.  In quality limits for instrument and attitude 
c    5.  In limits of acceptable date
c    6.  Not duplicates of science records (same time tags)
c    7.  Max dihedral temp not exceeded
c    8.  Accept sun, moon, earth limb angles within Min/Max limits.
c
c  throw bad records out of input buffer
c
c  if (total of science and sky records exceed maximum) then
c      signal error
c      set status to error 
c  endif
c
c  if (closeold and status) then
c      put only groups smaller than minimum in sky_buff
c  elseif (status)
c      put all records into sky_buff
c  endif
c
c  if (status) then
c      write summary of accepted and rejected records to report and log file
c  endif
c
c  fss_read_skymaps = status
c
c  return  
c  end (pdl)
c----------------------------------------------------------------------
C
C Changes:
C
C  D. Bouler, 12/19/91, SPR 9364 -- Increase number of eng records that
C  	can be read.
C
C  T. Hewagama, Add exclusion of data based on Sun, Moon, and Earth limb
C	angles. SER 9897 & 9898
C
C  T. Hewagama, 31-July-1992.  PIXLIST qualifier added to FSS.CLD, and
C	PLIST read terminator statement forced to quit at 100. SPR 9894
C
C  SPR 9942 - Add GALACTIC_LATITUDE to FSS output record, and select input 
C	skymap based on ABS( GALACTIC_LATITUDE ) > EXC_GALACTIC_LAT.
C	Tilak Hewagama, Hughes STX, September 1992
C----------------------------------------------------------------------

      implicit none
c
c Include files
c
      include       '(csa_pixel_input_rec)'
      include       '(csa_pixel_output_rec)'
      include       '(fss_include)'
      include       '(FUT_Params)'
c
c Function declarations
c
      integer*4      csa_read_pixels
      logical*4      time_gt
      real*8         aut_adt2dfloat
C
C Variable declarations
C
      integer*4      sky_recs              !Total skymap  recs read
      integer*4      sci_recs              !Total science recs read
      integer*4      i,j,k,l,m             !loop counters
      integer*4      status                !General return status
      integer*4      lu_input              !Input  skymap LU
      integer*4      num_in                !Number of input pixels
      integer*4      num_out               !Number of output pixels
      integer*4      max_recs              !Max number of pixel recs
      integer*4      blocks                !Block count for csa_read_pixels
      integer*4      recs                  !Number of skymap records
      integer*4      closeold              !Close out old skymap ?
      integer*4      minimum               !Minimum number for coadd group
      integer*4      num_bins              !Number of temperature bins
      integer*4      plist(*)              !List of pixels to process
      integer*4      mpmodes(*)            !List of modes to process
      integer*4      pixel                 !Current pixel number
      integer*4      mode                  !Current mode
      integer*4      inst_qual             !Instrument quality threshold
      integer*4      attit_qual            !Attitude quality threshold
      integer*4      lu                    !Generic LU
      integer*4      lu_rep                !Report file LU
      integer*4      chan                  !Channel being processed
      integer*4      start_chan            !Start channel
      integer*4      stop_chan             !Stop channel
      integer*4      chan_int              !Channel interval
      integer*4      out_bin               !Number rej for not in temp bin
      integer*4      out_mlist             !Number rej for not in mode list
      integer*4      out_mtmlist           !Number rej for not in mtm  list
      integer*4      out_plist             !Number rej for not in pixel list
      integer*4      out_qual              !Number rej for bad quality
      integer*4      out_time              !Number rej for time too early
      integer*4      out_dup               !Number rej for duplicate time 
      integer*4      out_closeold          !Number rej because of CLOSEOLD
      integer*4      out_dihed             !Number rej for high dihedral temp
      integer*4      num_pixels            !Number of pixels in list
      integer*4      num_modes             !Number of modes in list
      integer*4      time(2)               !Time (binary)
      integer*4      all_sky               !Total sky recs read in
      integer*4      earliest(2)           !Earliest acceptable time
      integer*4      chk_recs              !Number of records passing checks
      integer*4      groups                !Number of groups in input buffer
      integer*4      grp_size              !Size of current group
      integer*4      bptr                  !Temp bin number
      integer*4      length                !MTM scan length
      integer*4      speed                 !MTM scan speed
      integer*4      grp_mode              !Group mode
      integer*4      grp_pixel             !Group pixel
      integer*4      chknums(fss_max_buff) !Element numbers of good records
      integer*4      start_sky             !Start element number for sky recs
      integer*4      adt(2)                !Binary format time
      integer*4      inc                   !Increment between array elements in
                                           !dsci_times
      integer*4      time_index            !Science time array element that 
                                           !matches the skymap time (if any)
      integer*4      sci_times(2,4,fss_max_buff) ! Science times
      integer*4      sci_index             !Skymap index of recs from science
      integer*4      sky_index             !Skymap index of recs from skymap
      integer*4      stop_lu               !Max LU to print stats to
      integer*4      lu_int                !Interval for LUs
      integer*4      report                !Was report rewquested?
      integer*4      mtmlist(4)            !MTM modes allowed

      integer*4      out_sun               ! Recs outside sun angle range
      integer*4      out_moon              ! Recs outside moon angle range
      integer*4      out_earth             ! Recs outside earth limb angle range
      integer*4      out_galaxy            ! Recs outside galactic lat exclusion

      real*4         grp_temp              !Group ical temp
      real*4         bins(*)               !Temp bin midpoints
      real*4         temp                  !ical temp
      real*4         temptol               !Temperature tolerance
      real*4         max_dihed             !Max dihedral temp
      real*4 	     minsun                ! Minimum Sun angle
      real*4         maxsun                ! Maximum Sun angle
      real*4         minmoon               ! Minimum Moon angle
      real*4         maxmoon               ! Maximum Moon angle
      real*4         minearth              ! Minimum Earth limb angle
      real*4         maxearth              ! Maximum Earth limb angle
      real*4         sun                   ! Sun angle from record (deg)
      real*4         moon                  ! Moon angle from record (deg)
      real*4         earth                 ! Earth limb angle from record (deg)
      real*4         latitude              ! Galactic latitude
      real*4         Galaxy                ! Exclude recs with latitude < Galaxy
      
      real*8         dsci_times(fss_max_buff)  ! Science times
      real*8         dsky_time                 ! Skymap  time

      character*1    sci_pixdef            !Science pixel definition
      character*1    sky_pixdef            !Skymap  pixel definition      

      character*14   oldest                !Earliest acceptable time
      character*14   gmt                   !GMT format time
      character*14   early_gmt(100)        !Time tags earlier than OLDEST

      byte           level        /5/      !Index level for FIRAS

      logical*4      new_group             !Does current record start new group?
      logical*4      check_pixels          !Is pixel list in effect?
      logical*4      check_modes           !Is mode  list in effect?
      logical*4      in_bin                !Is temp in bin?
      logical*4      in_mlist              !Is mode in list?
      logical*4      in_mtmlist            !Is MTM mode in list?
      logical*4      in_plist              !Is pixel in list?
      logical*4      in_qual               !Is quality in bounds?
      logical*4      in_time               !Is time in bounds?
      logical*4      unique_time           !Time duplicated in science records?
      logical*4      pixdefs_match         !Does pixel definition match?
      logical*4      indices_match         !Do skymap indices match?

      logical*1      incsun                !Include record with this sun angle
      logical*1      incmoon               !Include record with this moon angle
      logical*1      incearth              !Include record with this earth angle
      logical*1      incgalaxy             !Include recs with this galactic lat.

c
c Declare externals
c
      external       fss_normal
      external       fss_allskyrej
      external       fss_duptime
      external       dsrch
      external       fss_pixelization
      external       fss_maxdata
c
c Dictionary and record declarations
c                                             
      dictionary 'fss_sssky'

      record /fss_sssky/ inp_buff(fss_max_buff)
      record /fss_sssky/ sky_buff(fss_max_buff)

      record /pixel_input_list/  inlist (6144)
      record /pixel_output_list/ outlist(6144)
c*
c******  Begin code *******************************************
c*
      fss_read_skymap = %loc(fss_normal)         
      status          = %loc(fss_normal)

      sky_recs     = 0
      chk_recs     = 0
      blocks       = 0
      out_bin      = 0
      out_mlist    = 0      
      out_plist    = 0     
      out_qual     = 0     
      out_time     = 0      
      out_closeold = 0
      out_dihed    = 0
      out_dup      = 0
      out_mtmlist  = 0
      all_sky      = 0
      start_sky    = sci_recs
      max_recs     = fss_max_buff
      num_in       = 6144
      inc          = 1

      out_sun = 0
      out_moon = 0
      out_earth = 0
      out_galaxy = 0

      call ct_gmt_to_binary(oldest, earliest)
c
c Convert science times to R*8 format so array can be searched
c for duplicate time tags from skymap using IMSL routine DSRCH
c
      do i = 1,sci_recs
         dsci_times(i) = aut_adt2dfloat(sci_times(1,chan,i))
      enddo
c
c Check to see if pixel list has been specified.
c
      if (plist(1) .ne. 9999) then
          num_pixels = 0
          do while (plist(num_pixels+1) .ne. -9999 .and. num_pixels.lt.100 )  
             num_pixels = num_pixels + 1
          enddo
          check_pixels  = .true.          
      else
          in_plist      = .true.
          check_pixels  = .false.          
      endif
c
c Check to see if microprocessor mode list has been specified.
c
      if (mpmodes(1) .ne. 9999) then
          num_modes = 0
          do while (mpmodes(num_modes+1) .ne. -9999)  
             num_modes = num_modes + 1
          enddo
          check_modes  = .true.
      else
          in_mlist     = .true.
          check_modes  = .false.
      endif
c
c Read in skymap to a temporary buffer and process.
c
      do i = 1, 6144
         inlist(i).pixel_no = i - 1 
         inlist(i).level_no = level
      enddo

      status = csa_read_pixels(lu_input, inlist,  num_in,  inp_buff,
     .                         max_recs, outlist, num_out, blocks)
      if (.not. status) call lib$signal(%val(status))
c
c Get total number of records read in.
c
      if (status) then
          recs = 0
          do i = 1, 6144
             recs = recs + outlist(i).no_records
          enddo
          all_sky = recs
      endif
c
c Check to see that pixelization is consistent.
c Must check the skymap index level also.
c
       if (status .and. (recs .gt. 0)) then

           sci_pixdef     = sky_buff(1).pixel_definition 
           sky_pixdef     = inp_buff(1).pixel_definition 
           pixdefs_match  = (sci_pixdef .eq. sky_pixdef)
           sci_index      = sky_buff(1).skymap_index
           sky_index      = inp_buff(1).skymap_index
           indices_match  = (sci_index .eq. sky_index)

           if ((.not. pixdefs_match) .or. (.not. indices_match)) then
                call lib$signal(fss_pixelization)
                status = %loc(fss_pixelization)
           endif
       endif
c
c Check records for agreement with thresholds.
c If record meets all requirements, increment
c the number of records passing the checks, 
c and save the array element number. 
c
      if (status .and. (recs .gt. 0)) then

          do i = 1,recs
c
c Check for pixel number 
c
             if (check_pixels) then 
                 in_plist = .false.
                 pixel = inp_buff(i).pixel_no
                 do j = 1,num_pixels
                    if (pixel .eq. plist(j)) in_plist = .true.
                 enddo
             endif
c
c Check for microprocessor mode 
c
             if (check_modes) then 
                 in_mlist = .false.
                 mode = inp_buff(i).sci_mode
                 do j = 1,num_modes
                    if (mode .eq. mpmodes(j)) in_mlist = .true.
                 enddo
             endif
c
c Check for MTM mode 
c
             speed  = inp_buff(i).mtm_scan_speed
             length = inp_buff(i).mtm_scan_length
             mode   = (2 * length) + speed
             in_mtmlist = .false.
             do j = 1,4
                if ((mode .eq. j-1) .and. (mtmlist(j) .eq. 1)) in_mtmlist=.true.
             enddo
c
c Check to see if the ical temp falls into a temp bin.
c
             temp   = inp_buff(i).ical_temp
             in_bin = .false.
             do j = 1,num_bins
                if ( abs(temp - bins(j)) .le. temptol) in_bin = .true.
             enddo
c
c Check to see if date of data is later than earliest allowable time
c
             if (time_gt(inp_buff(i).time, earliest)) then
                 in_time = .true.
             else
                 in_time = .false.
             endif
c
c Check to see if quality flags are within limits
c
             if ((inp_buff(i).data_quality     .le. inst_qual)    .and.
     &           (inp_buff(i).attitude_quality .le. attit_qual) ) then
                  in_qual = .true. 
             else
                  in_qual = .false.
             endif
c
c Check to see if time tag is duplicated in the science records.
c
             dsky_time  = aut_adt2dfloat(inp_buff(i).time)
             call dsrch(sci_recs, dsky_time, dsci_times, inc, time_index)
             if (time_index .gt. 0) then
                 call aut_dfloat2adt(dsky_time, adt)
                 call ct_binary_to_gmt(adt, gmt)
                 call lib$signal(fss_duptime, %val(1), gmt)
                 unique_time = .false.
             else
                 unique_time = .true.
             endif

c
c Check if sun, moon and earth limb angles are within Min/Max limits
c
             sun = Float( inp_buff(i).sun_angle ) * fac_att_conv
             moon = Float( inp_buff(i).moon_angle ) * fac_att_conv
	     earth = Float( inp_buff(i).earth_limb ) * fac_att_conv
	     latitude = Float( inp_buff(i).galactic_latitude ) * fac_att_conv

             if ( (sun .lt. minsun) .or. (sun .ge. maxsun) ) then
		incsun = .false.
	     else
             	incsun = .true.
             endif
             if ( (moon .lt. minmoon) .or. (moon .ge. maxmoon) )then
            	incmoon = .false.
	     else
             	incmoon = .true.
             endif
             if ( (earth .lt. minearth) .or. (earth .ge. maxearth) ) then
            	incearth = .false.
	     else
             	incearth = .true.
             endif

             if ( abs(latitude) .lt. Galaxy ) then
	    	incgalaxy = .false.
	     else
		incgalaxy = .true.
             endif


c
c Set counters for reasons to reject records.
c
             if (.not. in_qual) then
                  out_qual = out_qual + 1  
             elseif (.not. in_bin) then
                  out_bin  = out_bin  + 1
             elseif (inp_buff(i).dihedral_temp .gt. max_dihed) then
                  out_dihed = out_dihed + 1
             elseif (.not. in_plist) then
                  out_plist  = out_plist  + 1
             elseif (.not. in_mlist) then
                  out_mlist = out_mlist + 1
             elseif (.not. in_time) then
                  out_time = out_time + 1
                  if (out_time .le. 100) then
                     call ct_binary_to_gmt(inp_buff(i).time,early_gmt(out_time))
                  endif
             elseif (.not. unique_time) then
                  out_dup = out_dup + 1
             elseif (.not. in_mtmlist) then
                  out_mtmlist = out_mtmlist + 1
             elseif (.not. incsun) then
		out_sun = out_sun + 1
             elseif (.not. incmoon) then
            	out_moon = out_moon + 1
             elseif (.not. incearth) then
            	out_earth = out_earth + 1
             elseif (.not. incgalaxy) then
	    	out_galaxy = out_galaxy + 1
             endif

             if (in_bin  .and. in_mlist .and. in_plist        .and. 
     .           in_qual .and. in_time .and. unique_time      .and.
     .           (inp_buff(i).dihedral_temp .le. max_dihed)   .and.
     .           in_mtmlist .and. incsun .and. incmoon        .and.
     .		 incearth .and. incgalaxy)   then
 
                  chk_recs = chk_recs + 1
                  chknums(chk_recs) = i 

             endif

          enddo    !  (i = 1, recs)
      endif        !  (status .and. (recs .gt. 0)
c
c Replace the input array records with the records passing the checks.
c
      if ((chk_recs .gt. 0) .and. (chk_recs .lt. recs)) then 
           j = 0
           do i = 1, chk_recs
              j = j + 1                  
              inp_buff(j) = inp_buff(chknums(i))
           enddo
      endif
c
c Check that the total of science records and passed sky records 
c are less than the allowable maximum.
c
      if ( (sci_recs + chk_recs) .gt. fss_max_buff) then
            call lib$signal(fss_maxdata)
            status = %loc(fss_maxdata)
      endif
c
c If CLOSEOLD is specified, return only groups 
c less than the minimum coadd size in sky_buff.
c
      if (status .and. (.not. closeold) .and. (chk_recs .gt. 0)) then

          do i = 1, chk_recs
             sky_recs  = sky_recs  + 1
             start_sky = start_sky + 1
             sky_buff(start_sky) = inp_buff(i)
          enddo

      elseif (status .and. closeold .and. (chk_recs .gt. 0)) then
c
c Set parameters for current group equal to the first record.
c Use mode = (2 * length) + (speed); values are 0,1,2,3 (SS,SF,LS,LF).
c      
           bptr   = 1
           temp   = inp_buff(1).ical_temp
           in_bin = (abs(temp - bins(bptr)) .le. temptol)

           do while ( (.not. in_bin) .and. (bptr .lt. num_bins))
              bptr = bptr + 1
              in_bin = (abs(temp - bins(bptr)) .le. temptol)
           enddo

           grp_pixel  = inp_buff(1).pixel_no
           length     = 2 * inp_buff(1).mtm_scan_length  
           speed      = inp_buff(1).mtm_scan_speed
           grp_mode   = length + speed
           grp_temp   = bins(bptr) + temptol
           grp_size   = 0
c
c Group break will occur when:
c (1) Pixel number changes, (2) Scan mode changes, (3) temp crosses boundary.
c
           do i = 1, chk_recs
          
              length  = 2 * inp_buff(i).mtm_scan_length
              speed   = inp_buff(i).mtm_scan_speed 
              mode    = speed + length
              temp    = inp_buff(i).ical_temp
              pixel   = inp_buff(i).pixel_no

              new_group = (pixel .ne. grp_pixel) .or. 
     .                    (mode  .ne. grp_mode)  .or.
     .                    (temp  .gt. grp_temp)  

              if (.not. new_group) then
                   grp_size = grp_size + 1
              endif
c
c We found the beginning of a new group or reached the last record.
c If we are at the last record and it starts a new group, add the last
c record to the output buffer too.
c
               if (new_group) then
           
                   if (grp_size .lt. minimum) then
                       k = i - grp_size
                       l = i - 1
                       do j = k, l
                          sky_recs  = sky_recs  + 1
                          start_sky = start_sky + 1
                          sky_buff(start_sky) = inp_buff(j)
                       enddo
                   else
                       out_closeold = out_closeold + grp_size
                   endif
c
c If we are at the last record, add it to the output buffer.
c
                   if (i .eq. chk_recs) then
                       sky_recs  = sky_recs  + 1
                       start_sky = start_sky + 1
                       sky_buff(start_sky) = inp_buff(chk_recs)
                   endif

              elseif ((.not. new_group) .and. (i .eq. chk_recs)) then

                  if (grp_size .lt. minimum) then
                      k = i - grp_size + 1
                      do j = k, i
                         sky_recs  = sky_recs  + 1
                         start_sky = start_sky + 1
                         sky_buff(start_sky) = inp_buff(j)
                      enddo
                  else
                      out_closeold = out_closeold + grp_size
                  endif

              endif   !  (new_group .or. (i .eq. chk_recs))
c
c  Find the temperature bin to start on and reset the group parameters.
c
              if (new_group .and. (i .ne. chk_recs)) then 

                   bptr   = 1
                   in_bin = (abs(temp - bins(bptr)) .le. temptol)

                   do while ( (.not. in_bin) .and. (bptr .lt. num_bins))
                      bptr = bptr + 1
                      in_bin = (abs(temp - bins(bptr)) .le. temptol)
                   enddo

                   grp_size  = 1
                   grp_pixel = inp_buff(i).pixel_no
                   length    = 2 * inp_buff(i).mtm_scan_length 
                   speed     = inp_buff(i).mtm_scan_speed
                   grp_mode  = length + speed
                   grp_temp  = bins(bptr) + temptol

               endif    !   (new_group .and. (i .ne. chk_recs)) 
          enddo         !   while (i .le. chk_recs)
      endif             !   (status .and. .not. closeold .and. (recs .gt. 0))
c
c Reset transition flags to zero, as array will be resorted.
c
      if (status) then
          do i = sci_recs + 1, sci_recs + sky_recs
             sky_buff(i).transition_flag = 0
          enddo
      endif
c
c Write summary of records read and rejected to report file and log file.
c
      if (report) then
          stop_lu = lu_rep
          lu_int  = lu_rep - 6
      else
          stop_lu = 6
          lu_int  = 1
      endif

      do lu = 6, stop_lu, lu_int
         write(lu,fmt='(x,/,x,a,i1)') 
     .      'Input skymap rec failures for channel ',chan
         write(lu,fmt='(x,2a)') 'TOTSKY ACCEPT BADTIM PXLIST MDLIST TMPBIN ',
     .                          'MTMLST DATAQL CLSOLD DUPTIM DIHTMP'
         write(lu,10) all_sky, sky_recs, out_time, out_plist, out_mlist,
     .    out_bin, out_mtmlist,out_qual, out_closeold, out_dup, out_dihed
10       format(x,11(i6,x),/,x)
      enddo
c
c .................and the sun, moon, earth and galaxy faliures
      do lu = 6, stop_lu, lu_int
         write(lu,fmt='(x,a)') '   SUN   MOON  EARTH GALAXY'
         write(lu,11) out_sun, out_moon, out_earth, out_galaxy
11       format(x,4(i6,x),/,x)
      enddo
c
c Write out first 100 time tags that were earlier than OLDEST.
c
      if (out_time .gt. 0) then
          write(6,*) 'Time tags earlier than OLDEST:'
          i = min(100, out_time)
          do j = 1, i, 6
             k = min(i - j, 5)
             write(6,fmt='(x,6(a11,x))') (early_gmt(l)(:11),l=j,j+k)
          enddo
          write(6,*) ' '
       endif
c
c Signal error if all sky records were rejected.
c
      if (sky_recs .eq. 0) then
          call lib$signal(fss_allskyrej)
          status = %loc(fss_allskyrej)
      endif

      fss_read_skymap = status

      return
      end
