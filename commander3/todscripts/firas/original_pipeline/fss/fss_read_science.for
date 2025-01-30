      integer*4 function fss_read_science(sky_buff, sci_recs, galaxy, 
     .                                    pix_type, theta, phi, big_pixel,
     .                                    jstart, jstop, chan, plist, num_eng,
     .                                    sci_times, period, xcal,
     .                                    ical, refh, skyh, dihd, bolo, lu_rep,
     .                                    start_chan, stop_chan, chan_int,
     .                                    fdq_ext, report, minsun, maxsun,
     .                                    minmoon, maxmoon, minearth, maxearth,
     .                                    bins, num_bins, dbins, num_dbins, temptol)

C------------------------------------------------------------------------
C    PURPOSE: Reads the raw science data into an array.
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION: status = fss_read_science(sky_buff, sci_recs, galaxy, 
C    .                                     pix_type, theta, phi, big_pixel,
C    .                                     jstart, jstop, chan, plist,
C    .                                     num_eng, sci_times, period, xcal,
C    .                                     ical,refh, skyh, dihd, bolo, lu_rep, 
C    .                                     start_chan, stop_chan, chan_int,
C    .                                     fdq_ext, report, minsun, maxsun,
C    .                                     minmoon, maxmoon, minearth, maxearth,
C    .                                     bins, num_bins, dbins, num_dbins, temptol)
C
C    INPUT PARAMETERS:    R*4   galaxy                 Galactic plane exclusion
C                         I*4   pix_type               Pixelization chosen
C                         i*4   theta                  For scan_angle binning
C                         r*4   phi                    For orbit   binning
C                         i*4   big_pixel              For regular binning
C                         c*14  jstart                 Start time of processing
C                         c*14  jstop                  End   time of processing
C                         i*4   chan                   Current channel
C                         I*4   plist(*)               List of pixels to accept
C                         I*4   sci_times(2,4,max_buff) Science times
C                         i*4   num_eng(4)             Num eng recs/channel
C                         r*4   xcal(4,max_buff)       xcal temps
C                         r*4   ical(4,max_buff)       ical temps
C                         r*4   refh(4,max_buff)       ref horn temps
C                         r*4   skyh(4,max_buff)       sky horn temps
C                         r*4   dihd(4,max_buff)       dihedral temps
C                         r*4   bolo(4,max_buff)       bolometer temps
C                         I*4   period(2)              Orbital period (binary)
C                         I*4   lu_rep                 Report file LU
C                         I*4   start_chan             Start channel
C                         I*4   stop_chan              stop channel
C                         I*4   chan_int               channel interval
C                         I*4   report                 was report requested?
C			  R*4   minsun                 Minimum sun angle
C			  R*4   maxsun                 Maximum sun angle
C			  R*4   minmoon                Minimum moon angle
C			  R*4   maxmoon                Maximum moon angle
C			  R*4   minearth               Minimum earth angle
C			  R*4   maxearth               Maximum earth angle
C			  R*4   bins(100)              ICAL temperature bins
C			  R*4   dbins(100)             Dihedral temperature bins
C			  I*4   num_bins               Number of ICAL bins
C			  I*4   num_dbins              Number of dihedral bins
C			  R*4   temptol                ICAL temp. tolerance K.
C
C    OUTPUT PARAMETERS:   array sky_buff      array of short science records
C                         I*4   sci_rec       Number of science records read
C                         c*3   fdq_ext       FDQ file extent ('ED_' for edit)
C
C    SUBROUTINES CALLED:  None.
C
C    COMMON VARIABLES USED: None.
C 
C    INCLUDE FILES:   ct$library:ctuser.inc
C                     fut_params
C                     fss_include
C                     cct_filespec_fields_record      
C----------------------------------------------------------------------
C                     PDL FOR FSS_READ_SCIENCE
C                     D. BOULER, OCT, 1990
C
C  initialize variables
C
C  call lib$get_lun to get science logical unit number
C
C  adjust start and stop times by 30 minutes to avoid keyed read errors
C
C  open science archive
c
c  get FDQ file extension ('ED_' for edit)
c
c  do while ((status and (eng_ctr le num_eng(chan)) 
c
c     call ct_keyedread_arcv to use science time to get record
c
c     if (.not. status) then
c          signal error
c     else
c          increment total number of records read
c     endif
c
c     if (pixel list specified)   set pixel_in_list flag
c
c     if (science record is rejected) then
c         increment counter for reject reason
c     endif
c
c     if (status and pixel_in_list and (galactic lattitude in limits)) then
c
c         num_sci = num_sci + 1
c         
c         Fill in full science quauntities and eng quantities to short science
c
c         Set pixel number, pixel definition, and skymap index according
c         to pixelization scheme (regular, earth, scan mode, orbit)
c
c      endif
c
c      eng_ctr = eng_ctr + 1
c
c  enddo  !  { while status OK and counter < num_eng(chan) }
c
c  close full science archive
c
c  write summary of records read and rejected to report and log files
c
c  if (all science records rejected) then
c      signal error 
c      set status to error
c  endif
c
c  fss_read_science = status 
c
c  return
c  end (pdl)
C----------------------------------------------------------------------
C
C Changes:
C	Add exclusion of data based on Sun, Moon, and Earth limb angles.
C             Shirley M. Read, Hughes STX, February 1992. SERs 9897 & 9898
C
C	SPR 9540, FSS report file format is not clear to user
C		H. Wang, Hughes STX, March 4, 1992
C
C	PIXLIST qualifier added to FSS.CLD, and the PLIST read terminator
C		statement is forced to quit on reaching 100.
C		T. Hewagama, July 31, 1992, Hughes STX. SPR 9894
C
C	SPR 9942 - Add GALACTIC_LATITUDE to FSS output record, and select input 
C		skymap based on ABS( GALACTIC_LATITUDE ) > EXC_GALACTIC_LAT.
C		Tilak Hewagama, Hughes STX, September 1992
C
C	Store glitch count and galactic longitude in short science records
C               L. Rosen, Hughes STX, 10 February 1995
C
C       SPR 12197 - Modifications to implement /DBINS qualifier in place
C                   of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C----------------------------------------------------------------------

      implicit none
c
c Include files
c
      Include '($SSDef)'
      Include 'CT$LIBRARY:CTUser.Inc'
      Include '(FUT_Params)'
      Include '(fss_include)'
      Include '(cct_filespec_fields_record)'
c
c FUNCTION DECLARATIONS
c
      integer*4      lib$get_lun
      integer*4      sys$asctim
      integer*4      sys$bintim
      integer*4      lib$add_times
      integer*4      lib$sub_times
      integer*4      cct_get_filespec_fields

      logical*2      time_lt
      logical*2      time_gt
      logical*2      time_le
      logical*2      time_ge
c*
c* Declare externals
c*
      external       cct_get_filespec_fields
      external       ct_connect_read
      external       ct_binary_to_gmt
      external       ct_gmt_to_binary
      external       fss_badkeyread
      external       fss_nocatrecs
      external       fss_normal
      external       fss_badfsopen
      external       fss_badfsclose
      external       fss_noengrecs
      external       fss_allscirej
C
C Variable declarations
C
      integer*2	     ct_status(20)         ! COBETRIEVE status return
      integer*2      sample_mask           ! Mask to get start sample number
      parameter     (sample_mask = '0FFF'x)
      integer*2      gal_excl              !Galactic exclusion in 10-4 radians

      integer*4      start_sample          ! Start sample number
      integer*4      big_pixel             ! Used for regular pixelization
      integer*4      theta                 ! Used for scan_angle binning
      integer*4      pix_type              ! Type of pixelization chosen
      integer*4      start(2)              ! ADT start time
      integer*4      stop(2)               ! ADT stop time
      integer*4      plist(*)              ! List of pixels to accept
      integer*4      sci_times(2,4,fss_max_buff) ! full science times
      integer*4      Last_Orbit_Time(2)    ! Time of last science record 
      integer*4      Orbit_Gap(2)          ! Check orbit for gaps
      integer*4      period(2)             ! COBE orbital period
      integer*4      Time_Zero(2)/0,0/     ! zero time
      integer*4      Orbit_Pixel           ! Orbit oriented pixel number
      integer*4      i,j,k,l,ix            ! Counters
      integer*4      lu                    ! used to write to report & log file
      integer*4      status                ! General status return
      integer*4      lu_sci                ! LU for raw science
      integer*4      lu_eng                ! LU for engineering
      integer*4      ios                   ! I/O status return
      integer*4      adt(2)                ! Today's date in binary format
      integer*4      curr_time(2)          ! Current time in binary format
      integer*4      eng_ctr               ! Counter for eng recs 
      integer*4      num_eng(4)            ! Number of eng records/channel 
      integer*4      num_pixels            ! Number of pixels in list
      integer*4      sci_recs              ! Number of science records read
      integer*4      chan                  ! Current channel      
      integer*4      earliest(2)           ! Oldest useable date for data
      integer*4      bptr                  ! Pointer to current temp bin
      integer*4      time(2)               ! Time of short sci record
      integer*4      all_sci               ! total of science recs read
      integer*4      out_gal               ! number of recs out of gal exclusion
      integer*4      out_pix               ! number of recs not in pixel list
      integer*4      out_samp              ! Recs with start sample <> 1 
      integer*4      out_minsun            ! Recs with sun angle too low
      integer*4      out_maxsun            ! Recs with sun angle too high
      integer*4      out_minmoon           ! Recs with moon angle too low
      integer*4      out_maxmoon           ! Recs with moon angle too high
      integer*4      out_minearth          ! Recs with earth limb angle too low
      integer*4      out_maxearth          ! Recs with earth limb angle too high
      integer*4      pix_temp              ! temporary pixel number
      integer*4      lu_rep                ! Report file LU
      integer*4      start_chan            ! Start channel
      integer*4      stop_chan             ! Stop channel
      integer*4      chan_int              ! Channel interval
      integer*4      delta_adt(2)          ! Delta time
      integer*4      report                ! Was report requseted?
      integer*4      stop_lu               ! Max LU to write stats to
      integer*4      lu_int                ! Interval for LUs
   
      real*4         Latitude              ! Galactic latitude
      real*4         Scan_Angle            ! Scan angle in degrees
      real*4         Orbit_Phase           ! Orbit phase in degrees
      real*4         Last_Orbit_Phase      ! Orbit phase in degrees
      real*4         galaxy                ! Degrees to exclude near gal plane
      real*4         phi                   ! Used for orbit binning
      real*4         xcal(4,fss_max_buff)  ! xcal temps
      real*4         ical(4,fss_max_buff)  ! ical temps
      real*4         refh(4,fss_max_buff)  ! ref horn temps
      real*4         skyh(4,fss_max_buff)  ! sky horn temps
      real*4         dihd(4,fss_max_buff)  ! dihedral temps
      real*4         bolo(4,fss_max_buff)  ! bolometer temps
      real*4         temp                  ! Ical temp of record
      real*4 	     minsun                ! Minimum Sun angle
      real*4         maxsun                ! Maximum Sun angle
      real*4         minmoon               ! Minimum Moon angle
      real*4         maxmoon               ! Maximum Moon angle
      real*4         minearth              ! Minimum Earth limb angle
      real*4         maxearth              ! Maximum Earth limb angle
      real*4         sun                   ! Sun angle from record (deg)
      real*4         moon                  ! Moon angle from record (deg)
      real*4         earth                 ! Earth limb angle from record (deg)

      real*4         bins(100)             !ICAL temperature bins
      integer*4      num_bins              !Number of ICAL bins
      real*4         temptol               !ICAL temp. tolerance in K.
      real*4         dbins(100)            !Dihedral temperature bins
      integer*4      num_dbins             !Number of dihedral bins
      integer*2      bin_count         ! ICAL or dihedral bin index
      real*4         bin_max           ! ICAL or dihedral bin max temp
      integer*4      bin               ! ICAL bin of record
      integer*4      dbin              ! Dihedral bin of record
      logical*1      next              ! boolean to find ICAL or dihedral bin

      character*3    fdq_ext               ! FDQ file extension

      character*14   jstart                ! Start time (gmt)
      character*14   jstop                 ! Stop time  (gmt)
      character*14   sci_jstart            ! Science start time (gmt)
      character*14   sci_jstop             ! Science stop time  (gmt)
      character*14   key_time              ! Time used for keyed read
      character*14   gmt                   ! GMT time
      character*14   oldest                ! Oldest useable date for data

      character*16   delta_time            ! Time added to avoid keyed read err

      character*30   time_str              !Ascii current date and time

      character*128  sci_file              ! Name of full sci file

      logical*1      pixel_in_list         ! Is pixel in pixel list?
      logical*1      check_pixels          ! Are pixels being checked ?
      logical*1      incsun                !Include record with this sun angle
      logical*1      incmoon               !Include record with this moon angle
      logical*1      incearth              !Include record with this earth angle

c
c Declare records
c                                             
      dictionary 'fdq_sdf'
      dictionary 'fdq_eng'
      dictionary 'fss_sssky'

      record /filespec_fields/ fields_rec
      record /fss_sssky/       sky_buff(fss_max_buff) ! sci array
      record /fss_sssky/       init_sky               ! Short sci record
      record /fdq_sdf/         full                   ! Full sci records
c
c***  Begin code *******************************************
c
      fss_read_science = %loc(fss_normal)  
c
c Initialize variables
c
      status = %loc(fss_normal)

      if  (chan .ne. start_chan) then
           do i = 1, min(fss_max_buff, sci_recs)
              sky_buff(i) = init_sky
           enddo
      endif

      orbit_pixel = 0
      last_orbit_phase   = 999.0
      last_orbit_time(1) = 0
      last_orbit_time(2) = 0                   

      all_sci    = 0
      out_gal    = 0
      out_pix    = 0
      out_samp   = 0
      out_minsun = 0
      out_maxsun = 0
      out_minmoon = 0
      out_maxmoon = 0
      out_minearth = 0
      out_maxearth = 0

      gal_excl = iint(((galaxy * 10000.0 * 3.1415926) / 180.0) + 0.50)

      if ((plist(1) .ne. 9999) .and. (pix_type .le. 2)) then
          num_pixels = 0
          do while (plist(num_pixels+1) .ne. -9999 .and. num_pixels.lt.100 )  
             num_pixels = num_pixels + 1
          enddo
          check_pixels  = .true.          
      else
          check_pixels  = .false.          
          pixel_in_list = .true.
      endif
c
c The science times in the engineering records may overlap the start and stop
c times used to open the eng records. Subtract 30 minutes from JSTART and add
c 30 minutes to JSTOP to avoid getting keyed read errors when reading the 
c science records.
c
      delta_time = '   0 00:30:00.00'
      call sys$bintim(delta_time, delta_adt)

      call ct_gmt_to_binary(jstart, start)
      call ct_gmt_to_binary(jstop , stop)

      status = lib$add_times(stop, delta_adt, stop)
      if (.not. status) call lib$signal(%val(status))
      status = lib$sub_times(start, delta_adt, start)
      if (.not. status) call lib$signal(%val(status))

      call ct_binary_to_gmt(start, sci_jstart)
      call ct_binary_to_gmt(stop , sci_jstop )
c
c Open the full science archive.
c
      if (status) then 

          status = lib$get_lun (LU_sci)
          if (.not. status) call lib$signal(%val(status))

          sci_file = 'CSDR$FIRAS_RAW:FDQ_SDF_' // fac_channel_ids(chan) // 
     .               '/' // sci_jstart // ';' // sci_jstop // ';'

          if (status) then 
              open ( unit     = LU_sci,
     .               file     = sci_file,
     .               iostat   = ios,
     .               status   = 'OLD',
     .               useropen = CT_CONNECT_READ)
 
              if (ios .ne. 0) then
                  status = %loc(fss_badfsopen)
                  call lib$signal (fss_badfsopen, %val(1), %val(ios))
              end if
          endif
       endif     
c
c Read in the raw science. 
c The record will be used to build a short science entry in the sky_buff if
c     1. The lattitude is outside the galactic lattitudes being excluded 
c     2. The pixel is in the specified list.
c The type of pixelization, specified by pix_type, will determine the pixel.
c
      sci_recs = 0
      eng_ctr  = 1

      do while (status .and. (sci_recs .le. fss_max_buff) .and. 
     .         (eng_ctr .le. num_eng(chan)) )

         call ct_keyedread_arcv(, lu_sci, full, sci_times(1,chan,eng_ctr),
     .                            ct_status)

         if (ct_status(1) .ne. ctp_normal) then
             status = %loc(fss_badkeyread)    
             call ct_binary_to_gmt(sci_times(1,chan,num_eng(chan)), key_time)
             call lib$signal(fss_badkeyread, %val(3), %val(ct_status(1)),
     .                       %val(chan), key_time)
         else
             all_sci = all_sci + 1
         endif
c
c Check for pixel number if on command line.
c
         if (check_pixels .and. status .and. (pix_type .le. 2)) then 
             pixel_in_list = .false.
             if (pix_type .eq. 1) then
                 pix_temp = full.Attitude.Pixel_No/(4**(6-big_pixel))
             elseif (pix_type .eq. 2) then
                 pix_temp = full.Attitude.Terr_Pixel_No/(4**(6-big_pixel))
             endif
             do i = 1,num_pixels
                if (pix_temp .eq. plist(i)) pixel_in_list = .true.
             enddo
         endif

        Latitude     = full.Attitude.Galactic_Latitude * fac_att_conv
        start_sample = full.sci_head.sc_head16 .and. sample_mask
	sun = Float( full.attitude.sun_angle ) * fac_att_conv
        moon = Float( full.attitude.moon_angle ) * fac_att_conv
        earth = Float( full.attitude.earth_limb ) * fac_att_conv
c
c Set counters for reasons to reject science records.
c
         if (.not. pixel_in_list) then
               out_pix = out_pix + 1
         elseif (abs(latitude) .lt. galaxy) then
               out_gal  = out_gal  + 1
         elseif (start_sample .ne. 1) then
               out_samp = out_samp + 1
         endif
         incsun = .true.
         incmoon = .true.
         incearth = .true.
         if (sun .lt. minsun) then
            out_minsun = out_minsun + 1
            incsun = .false.
         endif
         if (sun .ge. maxsun) then
            out_maxsun = out_maxsun + 1
            incsun = .false.
         endif
         if (moon .lt. minmoon) then
            out_minmoon = out_minmoon + 1
            incmoon = .false.
         endif
         if (moon .ge. maxmoon) then
            out_maxmoon = out_maxmoon + 1
            incmoon = .false.
         endif
         if (earth .lt. minearth) then
            out_minearth = out_minearth + 1
            incearth = .false.
         endif
         if (earth .ge. maxearth) then
            out_maxearth = out_maxearth + 1
            incearth = .false.
         endif

        if (status .and. (abs(Latitude) .ge. Galaxy )  .and. 
     .      pixel_in_list .and. (start_sample .eq. 1)     
     .      .and. incsun .and. incmoon .and. incearth)  then

            sci_recs = sci_recs + 1
c
c Get the FDQ file extension.
c
            if (sci_recs .eq. 1) then
                status = cct_get_filespec_fields(lu_sci, fields_rec)
                if (.not. status) then
                     call lib$signal(%val(status))
                else
                     fdq_ext  = fields_rec.filename_extension(:3)
                endif
            endif

            sky_buff(sci_recs).Time(1)          = full.CT_Head.Time(1)
            sky_buff(sci_recs).Time(2)          = full.CT_Head.Time(2)
            sky_buff(sci_recs).Channel_Id       = full.Sci_Head.Chan_Id
            sky_buff(sci_recs).MTM_Scan_Speed   = full.Sci_Head.MTM_Speed
            sky_buff(sci_recs).MTM_Scan_Length  = full.Sci_Head.MTM_Length
            sky_buff(sci_recs).Adds_Per_Group   = full.Sci_Head.SC_HEAD9
            sky_buff(sci_recs).Sci_Mode         = full.Sci_Head.SC_HEAD1A
            sky_buff(sci_recs).data_quality     = full.dq_data.data_quality(110)
            sky_buff(sci_recs).attitude_quality = full.dq_data.data_quality(109)
            sky_buff(sci_recs).ICal_Temp        = ical(chan,eng_ctr)
            sky_buff(sci_recs).refhorn_temp     = refh(chan,eng_ctr)
            sky_buff(sci_recs).skyhorn_temp     = skyh(chan,eng_ctr)
            sky_buff(sci_recs).xcal_temp        = xcal(chan,eng_ctr)
            sky_buff(sci_recs).dihedral_temp    = dihd(chan,eng_ctr)
            sky_buff(sci_recs).bolometer_temp   = bolo(chan,eng_ctr)
            sky_buff(sci_recs).exc_galactic_lat = gal_excl
            sky_buff(sci_recs).sun_angle  = full.attitude.sun_angle
            sky_buff(sci_recs).moon_angle = full.attitude.moon_angle
            sky_buff(sci_recs).earth_limb = full.attitude.earth_limb
            sky_buff(sci_recs).galactic_latitude =
     .						full.attitude.galactic_latitude
            sky_buff(sci_recs).galactic_longitude =
     .                                         full.attitude.galactic_longitude
            sky_buff(sci_recs).glitch_count = full.sci_head.sc_head21
c
c Temporarily store 100 * ical bin number + dihedral bin number
c in the original_pixel field.
c After we sort by ical bin, dihedral bin, and time we then put
c original pixel number into this field.
c
            next = .true.
            bin = 0
            bin_count = 1
            do while (bin_count .le. num_bins .and. next)
               bin_max = bins(bin_count) + temptol
               if (sky_buff(sci_recs).ical_temp .le. bin_max) then
                  bin = bin_count
                  next = .false.                    ! ICAL bin found
               else
                  bin_count = bin_count + 1
               endif
            enddo
c
            next = .true.
            dbin = 0
            bin_count = 1
            do while (bin_count .le. num_dbins .and. next)
               bin_max = dbins(bin_count)
               if (sky_buff(sci_recs).dihedral_temp .le. bin_max) then
                  dbin = bin_count
                  next = .false.                    ! Dihedral bin found
               else
                  bin_count = bin_count + 1
               endif
            enddo
            sky_buff(sci_recs).original_pixel = 100 * bin + dbin
c
c Group data by sky based pixelization scheme
c
            if ( status .and. (pix_type .eq. 1) ) then
                 sky_buff(sci_recs).Pixel_No = 
     .           full.Attitude.Pixel_No/(4**(6-big_pixel))
                 sky_buff(sci_recs).Pixel_Definition = 'q'
                 sky_buff(sci_recs).Skymap_Index = Big_Pixel
c
c Group data by earth based pixelization scheme
c
             else if ( status .and. (pix_type .eq. 2) ) then
                  sky_buff(sci_recs).Pixel_No = 
     .            full.Attitude.Terr_Pixel_No/(4**(6-big_pixel))
                  sky_buff(sci_recs).Pixel_Definition = 'E'
                  sky_buff(sci_recs).Skymap_Index = Big_pixel
c
c Group data by scan angle
c
             else if ( status .and. (pix_type .eq. 3) ) then
                  Scan_Angle = full.Attitude.Scan_Angle * fac_att_conv
                  if (Scan_Angle .lt. 0.0 ) Scan_Angle = 360.0 + Scan_Angle
                  sky_buff(sci_recs).Pixel_No = int(Scan_Angle/theta)
                  sky_buff(sci_recs).Pixel_Definition = 'S'
                  sky_buff(sci_recs).Skymap_Index = Theta*10
c
c Group data by orbit phase angle
c
             else if ( status .and. (pix_type .eq. 4) ) then
                  Orbit_Phase = full.Attitude.Orbital_Phase*fac_att_conv
                  Orbit_Phase = Orbit_Phase - phi
                  if (Orbit_Phase .lt. 0.) Orbit_Phase=360. + Orbit_Phase
                  Call LIB$SubX(full.CT_Head.Time,Last_Orbit_Time,Orbit_Gap)
		  Call LIB$SubX(Time_Zero,Orbit_Gap,Orbit_Gap)
		  if ( (Orbit_Phase .lt. Last_Orbit_Phase) .or.
     .		        Time_LT(Orbit_Gap,period) )         then  
                         Orbit_Pixel = Orbit_Pixel + 1
                  end if
		  Last_Orbit_Phase = Orbit_Phase 
		  Last_Orbit_Time(1) = full.CT_Head.Time(1)
		  Last_Orbit_Time(2) = full.CT_Head.Time(2)
                  sky_buff(sci_recs).Pixel_No = Orbit_Pixel
                  sky_buff(sci_recs).Pixel_Definition = 'O'
                  sky_buff(sci_recs).Skymap_Index = Phi*10.

               end if      !   ( status OK and pix_type eq x )
           end if          !   ( XCAL out of horn and status OK and )
                           !   ( Quality is good and out of gal plane and )
                           !   ( Pixel in list and micro mode in list )
          
           eng_ctr = eng_ctr + 1

       enddo               !   ( while status OK and counter < MAX )
c
c Close the full science archive.
c
      if (status .ne. %loc(fss_badfsopen)) then
          call CT_CLOSE_ARCV (, LU_Sci, ct_status)
          if (ct_status(1) .ne. ctp_normal) then
              status = %loc(fss_badfsclose)
              call lib$signal(fss_badfsclose, %val(1), sci_file)
          endif
      endif
c
c Write summary of records read and rejected to report file and log file.
c
      if (status) then
          if (report) then
              stop_lu = lu_rep
              lu_int  = lu_rep - 6
          else
              stop_lu = 6
              lu_int  = 1
          endif

          do lu = 6, stop_lu, lu_int
             if (lu .eq. 6) then
                 write(lu,10) chan, all_sci,sci_recs
             else
                 write(lu,20) chan, all_sci,sci_recs
             endif          
             Write(lu,fmt='(x,a)') ' Records rejected for (sets may overlap):'
          
             write(lu,31)
          
             write(lu,30)  out_gal, out_pix, out_samp,
	1	     out_minsun, out_maxsun, out_minmoon,
	1	              out_maxmoon, out_minearth, out_maxearth 
10           format(x,/,x,' Channel = ',i1, 3x,'Total Science Records =',I7,
	1	3x, 'Records Accepted = ', I7/)
20           format('1',' Channel = ',i1, 3x,'Total Science Records =',I7,
	1	3x, 'Records Accepted = ', I7/)
31           format(21x,'Wrong',6x,'__(Less Than Min or Greater Than Max)__',/,
	1    1x,'Low Gal.   Not in   Start', 
	1	                6x,'___sun___      __moon___      __earth__',/,
	1    1x,'Latitude    Pixel  Sample',
	1	                6x,'min   max      min   max      min   max',/,
	1    1x,'--------   ------  ------',
	1	                6x,'----  -----    ----  -----    ----- -----')
30           format(x,i7,2x,i7,2x,i7,2x,i8,1x,i6,2x,i6,1x,i6,3x,i6,1x,i6)
     	       	
          enddo
       endif
c
c Signal error if all science records were rejected.
c
      if (sci_recs .eq. 0) then
          call lib$signal(fss_allscirej)
          status = %loc(fss_allscirej)
      endif

      fss_read_science = status

      return
      end
