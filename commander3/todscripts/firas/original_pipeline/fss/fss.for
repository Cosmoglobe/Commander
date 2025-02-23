      program fss

c*-----------------------------------------------------------------------------
c*
c*    Modified from fss to satisfy new requirements D. Bouler, STX, Sept, 1990.
c*
c                    PDL for FIRAS_SORT_SKY facility 
c                    D. Bouler, Jan, 1990 
c
c  Allocate virtual memory for arrays used to pass data
c
c  Call cut_display_banner and  cut_register_version to
c  establish which facility and which version is being run
c
c  Call fss_get_quals to get command line and qualifiers
c
c  if (status and RSE specified) call fss_get_rse to read RSE
c
c  If (report requested) then
c      Call fss_init_report to open report file 
c      Set fut_report_lun to report lun to send error messages to report
c  Endif
c
c  Call lib$establish to establish FUT_ERROR error handler
c
c  if (status) call fss_read_config to get grt switches and threshold values
c
c  if (status) then
c      call fss_read_eng to open eng archive, get science times and
c      controllable temps, and close eng archive.
c  endif
c
c  channel = start_channel
c  Do while ((channel .le. stop_channel) .and. status)
c 
c     if (status) then
c         call fss_read_science to open full science archive, fill
c         science and eng data into array of short science records, 
c         and close full science archive.
c     endif
c
c     if (status) then
c         call fss_open_skymaps to open the required skymaps, using either
c         default filenames or ones specified on command line
c     endif
c
c     if (status .and. (.not. FIRST_RUN)) then 
c         call fss_read_skymap to read in entire input skymap 
c     endif
c
c     if (status) then
c         call fss_sort to merge FDQ and sky data
c     endif
c
c     if (status) then
c         call fss_form_groups to form groups and set transition flags 
c     endif
c
c     if (status) then
c         call fss_write_skymaps to write new skymap 
c     endif
c
c     if (status) then
c         call fss_close_skymaps to close skymaps
c     endif
c
c     if (status and report requested) then
c         call fss_write_report to write report info for channel
c     endif
c
c     reset status to normal so processing can continue
c     increment channel
c  
c  enddo  !  (for each channel specified)
c
c  if (general return code is normal) then
c      signal normal completion of fes
c  else
c      signal fes abort
c  endif
c
c  close report file
c
c  end (pdl)
c
c
c
c    Modifications:
C 	Add exclusion of data based on Sun, Moon, and Earth limb angles.
C          SERs 9897 & 9898. Shirley M. Read, Hughes STX, February 1992
C
C	Add /PIXLIST option to pass in upto 100 pixels - SPR 9894.
C		Tilak Hewagama, Hughes STX, July 1992.
C
C	SPR 9942 - Add GALACTIC_LATITUDE to FSS output record, and select input 
C		skymap based on ABS( GALACTIC_LATITUDE ) > EXC_GALACTIC_LAT.
C				Tilak Hewagama, Hughes STX, September 1992
C
C Changes:
C	SPR 9943 - Report file includes CSDR$ archive references which
C			are not accessed by FSS.
C       February 1995 - Pass bins and num_bins to fss_write_skymap.  Pass
C            array of split flags from fss_form_groups to fss_write_skymap.
C
C       SPR 12197 - Modifications to implement /DBINS qualifier in place
C                   of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
c-----------------------------------------------------------------------------

      implicit none
c*
c* Include Files:
c*
      include 'CT$LIBRARY:CTUSER.INC' !Definitions for COBETRIEVE
      include '(FUT_PARAMS)'          !fac params
      include '(fss_include)'         !fss params
      include '($ssdef)'
      include '(fut_error)'
c*
c*    Externals:
c*
      external  fss_NORMAL              !Everything is OK
      external  fss_ABERR               !Everything is not OK
      external  FUT_Error
      external  lib$establish
c*
c*    Function Names:
c*
      integer*4 Cut_Register_Version
      integer*4 Cut_Display_Banner
      integer*4 fut_error               !Handles errors
      integer*4 lib$get_vm
      integer*4 lib$free_vm

      integer*4 fss_get_quals           !Get command line and qualifiers
      integer*4 fss_init_report         !Initializes the report
      integer*4 fss_get_rse             !Get the RSE if specified
      integer*4 fss_read_config         !Read in the config data
      integer*4 fss_read_eng            !Read in the eng data
      integer*4 fss_read_science        !Read in the new science data
      integer*4 fss_open_skymaps        !Open the skymaps
      integer*4 fss_read_skymap         !Reads the input skymap
      integer*4 fss_form_groups         !Form groups in ungrouped data
      integer*4 fss_sort                !Sort all data before write
      integer*4 fss_write_skymap        !Write the output skymap
      integer*4 fss_close_skymaps       !Close the skymaps
      integer*4 fss_write_report        !Writes the report
c*
c*    Declare Variables:
c*
      integer*4 Big_Pixel                  !User specified pixel resolution
      integer*4 Theta                      !User specified scan angle
      integer*4 pix_type                   !Type of pixelization
      integer*4 clen                       !Length of command line
      integer*4 num_eng(4)                 !Number of eng recs / channel
      integer*4 sky_recs                   !Number of recs from input skymap
      integer*4 closeold                   !Use stragglers from last month
      integer*4 write                      !Write to output archive or not
      integer*4 Instr_Qual                 !instrument data quality threshold
      integer*4 Attit_Qual                 !attitude data quality threshold
      integer*4 period(2)                  !COBE orbital period
      integer*4 tot_recs                   !Total number of recs , fdq and sky
      integer*4 Bad_ifgs                   !Total number of ifgs rejected
      integer*4 good_ifgs                  !Total number of good ifgs 
      integer*4 num_groups                 !Number of coadd groups formed
      integer*4 bad_groups                 !Number of coadd groups < minimum
      integer*4 sci_recs                   !Number of recs from full science 
      integer*4 chan                       !Channel 
      integer*4 start_chan                 !Channel to start processing
      integer*4 stop_chan                  !channel to stop  processing
      integer*4 chan_int                   !Interval between channels
      integer*4 First_Run                  !Flag indicating no input skymap
      integer*4 i,j,k,l                    !loop variables 
      integer*4 cnum                       !number of command line elements
      integer*4 ios                        !status of I/O operations
      integer*4 Status                     !general reteurn status 
      integer*4 qual_stat                  !Return from fss_get_quals
      integer*4 rep_stat                   !Return status for report
      integer*4 lu_config(2)               !LU for config files
      integer*4 LU_Output                  !LU for output file
      integer*4 lu_input                   !LU for input file
      integer*4 LU_Sci                     !LU for input archives
      integer*4 lu_rep                     !LU for report file
      integer*4 Minimum(4)                 !Minimum number of ifgs in a coadd 
      integer*4 plist(100)                 !List of pixels to be processed
      integer*4 mpmodes(10)                !Microprocessor modes to process
      integer*4 num_bins                   !Number of ICAL temp bins
      integer*4 num_dbins                  !Number of dihedral temp bins
      integer*4 report                     !Was report requested ?
      integer*4 vm_size                    !Amount of virtual memory (bytes)
      integer*4 mtmlist(4)                 !MTM modes (1=ss,2=sf,3=ls,4=lf)
                                           !(element = 1 if selected)
      integer*4 sky_jstart(2)              !JSTART of input skymap
      integer*4 fakeit                     !Allowable value for FAKEIT (0 or 1)
c
c Addresses for arrays in virtual memory.
c See call to get_vm for dimensions.
c     
      integer*4 sci_times(1)               !Sci times
      integer*4 pix_array (1)              !Pixels of coadd groups
      integer*4 mode_array(1)              !Modes of coadd groups
      integer*4 size_array(1)              !Sizes of coadd groups
      integer*4 bin_array (1)              !ICAL Temp bins of coadd groups
      integer*4 dbin_array (1)             !Dihedral Temp bins of coadd groups
      integer*4 split_array (1)            !Split flags of coadd groups
      integer*4 xcal(1)                    !xcal temps
      integer*4 ical(1)                    !ical temps
      integer*4 refh(1)                    !ref horn temps
      integer*4 skyh(1)                    !sky horn temps
      integer*4 dihd(1)                    !dihedral temps
      integer*4 bolo(1)                    !bolometer temps
      integer*4 sky_buff(1)                !Science and skymap recs

      real*4 bins(100)                     !ICAL Bin temps 
      real*4 dbins(100)                    !Dihedral Bin temps 
      real*4 Latitude                      !Galactic latitude
      Real*4 Phi                           !User specified orbit phase angle
      Real*4 Galaxy                        !User specified galactic latitude
      real*4 temptol                       !Max ICAL temp difference in a coadd 
      real*4 dtemptol                      !Max dihedral temp difference in a coadd 
      real*4 Scan_Angle                    !Scan angle in degrees
      real*4 Orbit_Phase                   !Orbit phase in degrees
      real*4 max_dihed                     !Max allowable temp for dihedrals
      real*4 minsun                        ! Minimum Sun angle
      real*4 maxsun                        ! Maximum Sun angle
      real*4 minmoon                       ! Minimum Moon angle
      real*4 maxmoon                       ! Maximum Moon angle
      real*4 minearth                      ! Minimum Earth limb angle
      real*4 maxearth                      ! Maximum Earth limb angle
 
      character*3  fdq_ext                 !FDQ file extension (e.g. 'ED_')

      character*10 detail                  !Detail level for report

      character*14 jstart                  !Start time for full science
      character*14 jstop                   !Stop  time for full science
      character*14 oldest                  !Oldest time for coaddition

      character*16 gmt                     !GMT time

      character*24 asc_orb                 !Ascii orbit time

      character*79 c_array(30)            !Array with command line

      character*128 sci_file               !Science file name 
      character*128 Repfile                !Report file name
      character*128 RSEfile                !RSE file name 
      character*128 inskymap               !Name of input skymap
      character*128 outskymap              !name of output skymap
      character*128 rse(16)     /16*' '/   !record selection expression
      character*128 engfile                !Name of eng input file

      logical*4     in_def                 !Was input skymap defaulted ?
      logical*4     out_def                !Was output skymap defaulted ?
      logical*4     eng_def                !Was eng file defaulted ?
      logical*4     sg_flag/.false./
c*
c* Declare record structures
c*
      dictionary 'fss_sssky'               !definition for FSS record
      dictionary 'fex_grtrawwt'            !definition for grt weights
      dictionary 'fex_grttrans'            !definition for grt trans temps

      record /fex_grtrawwt/  grtwt         ! grt weights
      record /fex_grttrans/  grttrans      ! grt transition temps

c*
c****   BEGIN CODE  **************************************************
c*
      status = %loc(fss_normal)

      call Cut_Register_Version(fss_version)
      call Cut_Display_Banner(6,80,'FIRAS Facility FSS_Sort_Sky')
      Write(6,fmt='(///)')

C*
C* GET VIRTUAL MEMORY FOR ARRAYS
C* SCIENCE TIMES  = (8 BYTES/TIME) * (4 CHANNELS) * (60,000 RECORDS/RUN)
C*
      vm_size = 8 * 4 * fss_max_buff
      status  = lib$get_vm( vm_size, sci_times(1)) 
      if (.not. status) then
           write(6,*) 'Insufficient VM for science time array' 
           call lib$show_vm()
      endif

C* QUANTITIES FOR COADD GROUPS = (4 BYTES) * (12,000 GROUPS/RUN)

      if (status) then
          vm_size =  4 * fss_max_groups
          status  =  lib$get_vm( vm_size, pix_array (1)) 
          status  =  lib$get_vm( vm_size, mode_array(1))
          status  =  lib$get_vm( vm_size, size_array(1))
          status  =  lib$get_vm( vm_size, bin_array (1))
          status  =  lib$get_vm( vm_size, dbin_array (1))
          vm_size =  2 * fss_max_groups
          status  =  lib$get_vm( vm_size, split_array (1))
          if (.not. status) then
               write(6,*) 'Insufficient VM for coadd arrays' 
               call lib$show_vm()
          endif
      endif

C* ENGINEERING QUANTITIES = (4 BYTES) * (4 CHANNELS) * (60,000 RECS/RUN)

      if (status) then
          vm_size =  4 * 4 * fss_max_buff
          status  =  lib$get_vm( vm_size, xcal(1)) 
          status  =  lib$get_vm( vm_size, ical(1)) 
          status  =  lib$get_vm( vm_size, refh(1))
          status  =  lib$get_vm( vm_size, skyh(1))
          status  =  lib$get_vm( vm_size, dihd(1))
          status  =  lib$get_vm( vm_size, bolo(1))
          if (.not. status) then
               write(6,*) 'Insufficient VM for eng arrays' 
               call lib$show_vm()
          endif
      endif

C* SKYMAP RECORDS = (skymap rec length) * (60,000 RECORDS / RUN)

      if (status) then
          vm_size =  fss_byte_recsize * fss_max_buff
          status  =  lib$get_vm( vm_size, sky_buff(1))
          if (.not. status) then
               write(6,*) 'Insufficient VM for sky rec array' 
               call lib$show_vm()
          endif
      endif
c*
c* GET COMMAND LINE AND QUALIFIERS
c*
      if (status) then
          qual_stat = FSS_get_quals (First_Run, start_chan, stop_chan,
     .                               chan_int, jstart, jstop,
     .                               Big_Pixel, Theta, Phi,
     .                               Galaxy, Instr_Qual, Attit_Qual, 
     .                               period, pix_type, report, cnum,
     .                               plist, mtmlist, mpmodes, oldest, 
     .                               closeold, write, c_array, inskymap,
     .                               outskymap, engfile, repfile, detail,
     .                               rsefile, bins, num_bins, dbins,
     .                               num_dbins, temptol, 
     .                               in_def, out_def, eng_def, fakeit,
     .                               minsun, maxsun, minmoon,
     .                               maxmoon, minearth, maxearth )
      endif
c*
c* Maximum Dihedral Temperature from DBINS
c*
      max_dihed = dbins(num_dbins)
c*
c* GET RSE
c*
      if (status .and. (rsefile .ne. ' ')) then
          status = FSS_get_rse ( rsefile, RSE )
          if (.not. status) call lib$signal(%val(status))
      endif
c*
c* Initialize report
c*
      if (report) then
          status = fss_init_report(repfile, lu_rep, c_array, qual_stat,
     .                             cnum, rse, rsefile)
          if (status) call lib$establish(fut_error)
      endif
c*
c* READ CONFIG DATA
c*
      if (status) then
          status = fss_read_config(minimum, grtwt, grttrans)
      endif
c*
c* READ IN ENGINEERING DATA AND GET TIMES FOR SCIENCE RECORDS
c*
      if (status) then 
          status = fss_read_eng( jstart, jstop, rse, num_eng, 
     .                           %val(xcal(1)), %val(ical(1)), %val(refh(1)),
     .                           %val(skyh(1)), %val(dihd(1)), %val(bolo(1)),
     .                           %val(sci_times(1)), grtwt, grttrans,  
     .                           start_chan, stop_chan, chan_int, lu_rep,
     .                           instr_qual, attit_qual, mpmodes,
     .                           mtmlist,oldest,bins,num_bins,temptol,
     .                           max_dihed, engfile, fakeit, report)
      endif

c*
c* DO FOR EACH CHANNEL
c*
      chan = start_chan

      do while ((chan .le. stop_chan) .and. status)
c*
c* READ IN ALL RAW SCIENCE DATA
c*
          if (status) then 
              status = fss_read_science(%val(sky_buff(1)),sci_recs, galaxy, 
     .                                  pix_type, theta, phi, big_pixel,
     .                                  jstart, jstop, chan, plist,
     .                                  num_eng, %val(sci_times(1)), period, 
     .                                  %val(xcal(1)), %val(ical(1)),
     .                                  %val(refh(1)), %val(skyh(1)), 
     .                                  %val(dihd(1)), %val(bolo(1)), lu_rep,
     .                                  start_chan, stop_chan, chan_int,
     .                                  fdq_ext, report, minsun, maxsun,
     .                                  minmoon, maxmoon, minearth, maxearth,
     .                                  bins, num_bins, dbins, num_dbins,
     .                                  temptol)

          endif
c*
c* OPEN SKYMAP FILES
c*
         if (status) then
             status = fss_open_skymaps(lu_input, lu_output, first_run,
     .                                 closeold, oldest, 
     .                                 inskymap, outskymap, write,
     .                                 jstart, jstop, chan, fdq_ext,
     .                                 sky_jstart, in_def, out_def, engfile)
         endif
c*
c* READ IN ENTIRE INPUT SKYMAP
c*
         if (status .and. (.not. first_run)) then
             status = fss_read_skymap(lu_input, %val(sky_buff(1)), 
     .                                sci_recs, sky_recs, closeold,
     .                                minimum(chan), bins, num_bins, dbins, 
     .                                num_dbins, mpmodes, plist, instr_qual,
     .                                attit_qual, oldest, lu_rep, chan, 
     .                                start_chan, stop_chan, chan_int, temptol,
     .                                %val(sci_times(1)), max_dihed ,report,
     .                                mtmlist, galaxy, minsun, maxsun,
     .                                minmoon, maxmoon, minearth, maxearth)

         endif
c*
c* MERGE AND SORT ALL DATA 
c*
         if (status) then
             status = fss_sort(%val(sky_buff(1)), sci_recs, 
     .                         sky_recs, tot_recs, chan)
         endif
c*
c* FORM COADD GROUPS FROM MERGED UNCOADDED DATA
c*
         if (status) then
             status = fss_form_groups(%val(sky_buff(1)), tot_recs, chan,
     .                                num_groups, bad_groups, good_ifgs, 
     .                                bad_ifgs, bins, num_bins, dbins, num_dbins,
     .                                temptol, minimum, %val(pix_array(1)),
     .                                %val(mode_array(1)), %val(bin_array(1)),
     .                                %val(dbin_array(1)), %val(size_array(1)),
     .                                %val(split_array(1)))
          endif
c*
c* WRITE OUTPUT SKYMAP
c*
         if (status .and. write) then
             status = fss_write_skymap (lu_output, %val(sky_buff(1)), tot_recs,
     .                            bins, num_bins, dbins, num_dbins, temptol, 
     .                            %val(split_array(1)), %val(size_array(1)))
         endif
c*
c* CLOSE OUTPUT SKYMAP
c*
         if (status) then
             status = fss_close_skymaps(lu_input, lu_output, first_run, write,
     .                                  closeold, jstart, jstop, sky_jstart,
     .                                  in_def, eng_def)
         endif
c*
c* WRITE THE SUMMARY REPORT INFO FOR THE CHANNEL
C*
         if (report .and. status) then 
             status = fss_write_report(num_groups, bad_groups, good_ifgs,
     .                                 bad_ifgs, %val(pix_array(1)),
     .                                 %val(mode_array(1)),%val(bin_array(1)), 
     .                                 %val(dbin_array(1)), %val(size_array(1)), 
     .                                 minimum, temptol, detail, chan, 
     .                                 lu_rep, inskymap, outskymap, 
     .                                 engfile, first_run, write)

         endif
c*
c* INCREMENT CHANNEL AND RESET STATUS
c*
         if (chan .lt. stop_chan) status = %loc(fss_normal)
         chan = chan + chan_int

       enddo      !   while ((chan .le. stop_chan) .and. status)
C*
C*  FREE VIRTUAL MEMORY
C*
      vm_size = 8 * 4 * fss_max_buff
      call lib$free_vm( vm_size, sci_times(1)) 

      vm_size = 4 * fss_max_groups
      call lib$free_vm( vm_size, pix_array (1)) 
      call lib$free_vm( vm_size, mode_array(1))
      call lib$free_vm( vm_size, size_array(1))
      call lib$free_vm( vm_size, bin_array (1))
      call lib$free_vm( vm_size, dbin_array (1))
      vm_size = 2 * fss_max_groups
      call lib$free_vm( vm_size, split_array(1))

      vm_size = 4 * 4 * fss_max_buff
      call lib$free_vm( vm_size, xcal(1)) 
      call lib$free_vm( vm_size, refh(1))
      call lib$free_vm( vm_size, skyh(1))
      call lib$free_vm( vm_size, dihd(1))
      call lib$free_vm( vm_size, bolo(1))

      vm_size = fss_byte_recsize * fss_max_buff
      call lib$free_vm( vm_size, sky_buff(1))

c*
c* SIGNAL COMPLETION STATUS
c*
       if (status) then
 	  call lib$signal(fss_normal)
       else
 	  call lib$signal(fss_aberr)
       endif
c*
c* CLOSE REPORT FILE
c*
      if (report) then
          close(lu_rep,iostat=ios)
          if (ios .ne. 0) then
              call errsns(ios,,,,status)
              call lib$signal(%val(status))
          endif
          fut_report_lun = 0
      endif

       if (status) then
          call exit(ss$_normal)
       else
          call exit(ss$_abort)
       endif

       end
