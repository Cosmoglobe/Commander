      integer*4 function FSS_get_quals (First_Run, start_chan, stop_chan,
     .                                  chan_int, jstart, jstop,
     .                                  Big_Pixel, Theta, Phi,
     .                                  Galaxy, Instr_Qual, Attit_Qual, 
     .                                  period, pix_type, report, cnum,
     .                                  plist, mtmlist, mpmodes, oldest, 
     .                                  closeold, write, c_array, inskymap,
     .                                  outskymap, engfile, repfile, detail,
     .                                  rsefile, bins, num_bins, dbins,
     .                                  num_dbins, temptol,
     .                                  in_def, out_def, eng_def, 
     .                                  fakeit, minsun, maxsun,
     .                                  minmoon, maxmoon, minearth, maxearth )

C------------------------------------------------------------------------
C
C    PURPOSE: This function gets the command line qualifiers, except for RSE.
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION: status = FSS_get_quals (First_Run, start_chan, stop_chan,
C     .                                  chan_int, jstart, jstop,
C     .                                  Big_Pixel, Theta, Phi,
C     .                                  Galaxy, Instr_Qual, Attit_Qual, 
C     .                                  period, pix_type, report, cnum,
C     .                                  plist, mtmlist, mpmodes, oldest, 
C     .                                  closeold, write, c_array, inskymap,
C     .                                  outskymap, engfile, repfile, detail,
C     .                                  rsefile, bins, num_bins, dbins, 
C     .                                  num_dbins, temptol,  
C     .                                  in_def, out_def, eng_def,
C     .                                  fakeit, minsun, maxsun,
C     .                                  minmoon, maxmoon, minearth, maxearth)
C  
C    INPUT PARAMETERS:     NONE.   
C
C    OUTPUT PARAMETERS:    i*4    first_run         Is this first run of FSS?
C                          i*4    start_chan        Start channel
C                          i*4    stop_chan         Stop  channel
C                          i*4    chan_int          Interval between channels
C                          c*14   jstart            start time
C                          c*14   jstop             stop time
C                          i*4    big_pixel         Value for big pixel
C                          i*4    theta             Value for scan angle
C                          r*4    phi               Orbit angle
C                          r*4    galaxy            Galactic exclusion angle
C                          i*4    instr_qual        Instrument quality
C                          i*4    attit_qual        Attitude quality
C                          i*4    period(2)         ADT format orbital period
C                          i*4    pix_type          Pixel type
C                          i*4    report            Was report requested ?
C                          i*4    plist(100)        List of pixels to process
C                          i*4    cnum              Num lines for command line
C                          i*4    mtmlist(4)        MTM modes to process
C                          i*4    mpmodes(10)       List of micro modes
C                          c*14   oldest            Earliest time for data
C                          I*4    closeold          Read in only uncoadded ifgs
C                          I*4    write             Write to archive
C                          c*79   c_array(30)       Command line
C                          c*128  inskymap          Input skymap name
C                          c*128  outskymap         Output skymap name
C                          c*128  engfile           Eng file name
C                          c*128  repfile           Report file name
C                          c*10   detail            Level of report detail
C                          C*128  rsefile           Rse file name
C                          R*4    bins(100)         ICAL Temp bins
C                          R*4    dbins(100)        Dihedral Temp bins
C                          I*4    num_bins          Number of ICAL temp bins
C                          I*4    num_dbins         Number of Dihedral temp bins
C                          R*4    temptol           Tolerance for ICAL temp bins
C                          L*4    in_def            Was input skymap defaulted?
C                          L*4    out_def           Was output skymap defaulted?
C                          I*4    fakeit            FAKEIT data allowed?
C			   R*4    minsun            Minimum sun angle
C			   R*4    maxsun            Maximum sun angle
C			   R*4    minmoon           Minimum moon angle
C			   R*4    maxmoon           Maximum moon angle
C			   R*4    minearth          Minimum earth angle
C			   R*4    maxearth          Maximum earth angle
C                    
C    SUBROUTINES CALLED:   cli$present
C                          upm_get_float
C                          upm_get_longword
C                          fut_orbital_period
C                          cli$get_value
C                          upm_get_value
C
C    COMMON VARIABLES USED:  None.
C 
C    INCLUDE FILES:          '($ssdef)'
C                            'CT$LIBRARY:CTUSER.INC'      
C
C----------------------------------------------------------------------
C
C Changes:
C
c    SER 7126  Rewrite of FSS to meet new requirements.
c              D. Bouler, STX, OCT, 1990
C
c    Add qualifiers for exclusion of data by Sun, Moon, and Earth limb angles.
c	       Shirley M. Read, Hughes STX, February 1992. SERs 9897 & 9898
c
c    SER 9538 FSS problem with MTM command line qulifier H. Wang, Hughes STX,
c		March, 1992
c
c    PIXLIST qualifier added, Tilak Hewagama, Hughes/STX, 31-July-1992. SPR 9894
C
C
C    SPR 12197 - Modifications to implement /DBINS qualifier in place
C                of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C----------------------------------------------------------------------
c
c                  PDL for FSS_GET_QUALS
c                  D. Bouler, Nov, 1990
c
c   Set return status to fss_normal
c
c   Call CT_INIT to set up COBETRIEVE
c   Get JSTART and JSTOP values and form time range 
c   Get current date and time
c
c   Get whether report was requested
c   if (report .and. report name specified on command line) then
c       get report file name
c   elseif (report .and. report name not specified on command line) then
c       build default report file name
c   endif
c   Get level of detail for report
c
c   Get the type of pixelization:
c   pix_type = 1 for regular pixelization, big_pixel = 6
c   pix_type = 2 for earth-based pixels.   big_pixel = 6
c   pix_type = 3 for scan-angle pixels,    theta     = specified angle
c   pix_type = 4 for orbital angle pixels, phi       = specified angle
c
c   Get value for galactic exclusion
c   Get instrument and attitude quality thresholds
c   Get orbital period if orbital pixels are used
c   Set flag for first_run option
c   Set flag for closeold option
c   Set flag for write option
c   Set flag for fakeit option
c   Get the channels to process   
c   Get MTM modes to process
c   Get pixel list to process
c   Get microprocessor mode list to process
c   Get time for oldest data to process (default = 89320)
c   Get input skymap file extension, and whether it was defaulted
c   Get output skymap file extension, and whether it was defaulted
c   Get eng file extension, and whether it was defaulted
c   Get RSE file name 
c   Get temperature tolerance
c   Get Max dihedral Temp
c   Get the minimum and maximum Sun, Moon, and Earth limb angle exclusions
c   Get ICAL Temperature bins
c   Get Dihedral Temperature bins
c   Check for dihedral bin sorting
c   Get command line for FSS run
c   return
c   end (pdl)
c
c-----------------------------------------------------------------------------

      implicit none
c*
c*    Include Files:
c*
      include '($ssdef)'
      include '(fut_params)'
      include 'CT$LIBRARY:CTUSER.INC'      ! COBETRIEVE definitions
c*
c*    Functions:
c*
      integer*4 fut_orbital_period
      integer*4 cli$present                !Present on command line?
      integer*4 upm_get_longword
      integer*4 upm_get_float
      integer*4 upm_get_value
      integer*4 cli$get_value
c*
c*    Externals:
c*
      external fss_ctiniterr               !CT-INIT failure
      external CLI$_Present                !qualifier present
      external CLI$_defaulted              !qualifier defaulted
      external CLI$_negated                !qualifier negated
      external CLI$_absent                 !qualifier absent
      external fss_normal                  !Everything OK signal
      external fss_bin_overlap             !ICAL Temp bins on command line overlap
      external fss_dihed_sort              !Dihed Temp bins on command line not sorted
      external upm_Comma    
      external upm_Concat   
      external upm_Pres     
      external upm_Negated  
      external upm_Locpres  
      external upm_Locneg   
      external upm_Defaulted
      external upm_Absent   
      external upm_Invnum   
c*
c*    Variables:
c*
      Parameter Absent_upm = '0B6DBE80'X   !Upm value absent

      integer*2 ct_status(20)              !COBETRIEVE return array

      integer*4 Big_Pixel                  !User specified pixel resolution
      integer*4 Theta                      !User specified scan angle
      integer*4 pix_type                   !Parameter for grouping the Ifgs
      integer*4 write                      !Write output records or not
      integer*4 closeold                   !Use stragglers from last month
      integer*4 i,j,k,l                    !Loop counters
      integer*4 slash(100)                 !Position of slashes in comm line
      integer*4 num_slash                  !Number of slashes in command line
      integer*4 line_len      /78/         !Max length of printed command line 
      integer*4 pix_line_len  /74/         !Max len of printed com line for pixs
      integer*4 use_line_len               ! value = (pix_line_len or line_len)
      integer*4 qual_len                   !Length of qualifier
      integer*4 line_start                 !Start of command line element
      integer*4 mpmodes(10)                !Microprocessor modes to process
      integer*4 ios/0/                     !Input/Output status
      integer*4 status                     !general status
      integer*4 rep_stat                   !Status return for report file name
      integer*4 length                     !number of characters in Value
      integer*4 First_Run                  !Flag for whether holding files exist
      Integer*4 period(2)                  !Orbital period in ADT.
      Integer*4 Time_Zero(2)/0,0/          !Convert to a delta time.
      integer*4 Instr_Qual                 !Instrument quality threshold
      integer*4 Attit_Qual                 !Attitude quality threshold
      integer*4 plist(100)                 !List of pixels to process
      integer*4 start_chan                 !Channel to start processing
      integer*4 stop_chan                  !Channel to stop  processing
      integer*4 chan_int                   !Interval between channels
      integer*4 c_len                      !Length of command line
      integer*4 adt(2)                     !Today's date (binary)
      integer*4 num_bins                   !Number of ICAL temp bins
      integer*4 num_dbins                  !Number of dihedral temp bins
      integer*4 num_mtm                    !Number of mtm modes
      integer*4 cnum                       !Number command line elements 
      integer*4 report                     !Was report requested ?
      integer*4 mtmlist(4)                 !List of MTM modes (1=process)
                                           !(1=ss,2=sf,3=ls,4=lf)
      integer*4 fakeit                     !Is FAKEIT data allowed?
      integer*4 rse_stat                   !Return from getting RSE name

      Real*4 Phi                           ! User specified orbit phase angle
      Real*4 Galaxy                        ! User specified galactic latitude
      Real*4 R4_Orbital_PD                 ! User specified orbital pd in min.
      real*4 bins(100)                     ! ICAL Temp bins
      real*4 dbins(100)                    ! Dihedral Temp bins
      real*4 temptol                       ! ICAL Temp tolerance      
      real*4 rtemp                         ! Temporary variable
      real*4 minsun                        ! Minimum Sun angle
      real*4 maxsun                        ! Maximum Sun angle
      real*4 minmoon                       ! Minimum Moon angle
      real*4 maxmoon                       ! Maximum Moon angle
      real*4 minearth                      ! Minimum Earth limb angle
      real*4 maxearth                      ! Maximum Earth limb angle

      Real*8 R8_Orbital_PD                 ! Orbital pd in 10^-7 sec.

      character*3   mtm(4)                 !MTM mode (SS,SF,LS,LF,ALL)

      character*4   asc_i1                 !Ascii version of integer            
      character*4   asc_i2                 !Ascii version of integer            
      
      character*5   channel                !Channel from command line
      character*8   asc_f                  !Ascii version of R*4 number

      character*6   atemptol               !Ascii version of ICAL temp tolerance
      
      character*6   aminsun                !Ascii version of min sun
      character*6   amaxsun                !Ascii version of max sun
      character*7   aminmoon               !Ascii version of min moon
      character*7   amaxmoon               !Ascii version of max moon
      character*8   aminearth              !Ascii version of min earth
      character*8   amaxearth              !Ascii version of max earth

      character*10  detail                 !Level of detail for report

      character*14  Oldest                 !Earliest time to process ifgs
      character*14  today                  !Today's date (Gmt)
      character*14  jstart                 !Start time (gmt)
      character*14  jstop                  !Stop  time (gmt)

      character*79  c_array(30)            !Holds command line

      character*128 inskymap               !Input skymap name
      character*128 outskymap              !Output skymap name
      character*128 engfile                !Eng file name
      character*128 repfile                !Report file name      
      character*128 rsefile                !RSE file name

      character*255 Value                  !Generic string for cli$get_value

      character*2048 c_line                !Command line for FSS invocation

      logical*4      repdef                !Was report name defaulted ?
      logical*4      bins_defaulted        !Were temp bins defaulted  ?
      logical*4      overlap               !Do ICAL temp bins overlap ?
      logical*4      dihed_sort            !Are Dihedral temp bins unsorted ?
      logical*4      in_def                !Was input  skymap defaulted ?
      logical*4      out_def               !Was output skymap defaulted ?
      logical*4      eng_def               !Was eng file defaulted ?

c* These definitions for pixlist qualifier ( T.H. 31-July-92 )
      integer ilis
      integer pixunit
      integer*4 lib$get_lun
      integer*4 iosp/0/                    !Input/Output status pixels read
      integer*4 statusp                    !General status for pixel to be read
      integer*4 plists(100)		   !Temporary storage for plist
      namelist /pixlist/ plists


c*
c***  BEGIN CODE  ******************************************************
c*
      fss_get_quals = %loc(fss_Normal)
      status        = %loc(fss_Normal)

c*
c* Get todays date and time
c*
      call sys$gettim(adt)
      call ct_binary_to_gmt(adt, today)
c*
c* CALL CT_Init TO PREPARE COBETRIEVE
c*
      call ct_init (ct_status)
      if (ct_status(1) .ne. ctp_normal) then
          call lib$signal(%val(ct_status(1)))
          status = %loc(fss_ctiniterr)
      endif
c*
c* Start command line; build as each qualifier is found
C*
      c_line = 'FSS'
      c_len  = 3
c*
c* Get jstart and jstop 
c*
      status = cli$get_value ('JSTART', jstart, length)
      if (.not. status) then
           call lib$signal(%val(status))
      else      
           if (length .lt. 14) then
               jstart(length+1:14) = fac_jstart_default(length+1:14)
           endif
           c_line = c_line(:c_len) // '/JSTART=' // jstart
           c_len  = c_len + 22 
      endif
 
      if (status) then
          status = cli$get_value ('JSTOP', jstop, length)
          if (.not. status) then
               call lib$signal(%val(status))
          else
               if (length .lt. 14) then 
                   jstop(length+1:14) = fac_jstop_default(length+1:14)
               endif
               c_line = c_line(:c_len) // '/JSTOP=' // jstop
               c_len  = c_len + 21 
          endif
       endif
c
c Get report qualifiers in command line array. 
c
      report = cli$present('REPORT')
      repdef = report .eq. %loc(cli$_defaulted)

      if (status .and. report .and. (.not. repdef)) then
          rep_stat = cli$get_value('REPORT', repfile, length)
          if ((.not. status) .and. (rep_stat .ne. %loc(cli$_absent))) then
              call lib$signal(%val(rep_stat))
          endif
          if (length .eq. 0) repdef = .true.
      endif       

      if (status .and. report .and. repdef) then
          repfile = 'FSS_'//jstart(:7)//'_'//jstop(:7)//'.REP_'//today(:9)
          call str$trim(repfile, repfile, length)
      endif

      if (status .and. (.not. report)) then
           c_line = c_line(:c_len) // '/NOREPORT'
           c_len  = c_len + 9
      elseif (status .and. report) then
          c_line = c_line(:c_len) // '/REPORT=' // repfile(:length)
          c_len  = c_len + 8 + length
      endif          
c
c Get Detail level for report
c
      if (status .and. report) then
          if (cli$present('DETAIL') .eq. %loc(cli$_defaulted)) then
              detail = 'SUMMARY'
              c_line = c_line(:c_len) // '/DETAIL=' // detail(:7)
              c_len  = c_len + 15
          else
              detail = 'FULL'
              c_line = c_line(:c_len) // '/DETAIL=' // detail(:4)
              c_len  = c_len + 12
          endif
      endif          
c*
c* Get the type of pixelization and values for scan angle and orbit
c*
      if (cli$present('TYPEPIX') .and. status) Then
          if (cli$present('SKY')) Then
             pix_type = 1
             status = upm_get_longword('SKY',big_pixel)
             write(ASC_I1,fmt='(i1)') big_pixel
             c_line = c_line(:c_len) // '/TYPEPIX=SKY=' // ASC_I1
             c_len  = c_len + 14
          else if (cli$present('EARTH')) Then
             pix_type = 2
             status = upm_get_longword('EARTH',big_pixel)
             write(ASC_I1,fmt='(i4)') big_pixel
             c_line = c_line(:c_len) // '/TYPEPIX=EARTH=' // ASC_I1
             c_len  = c_len + 19
          else if (cli$present('SCAN_ANGLE')) Then
             pix_type = 3
             status = upm_get_longword('SCAN_ANGLE',theta)
             write(ASC_I1,fmt='(i4)') theta
             c_line = c_line(:c_len) // '/TYPEPIX=SCAN=' // ASC_I1
             c_len  = c_len + 18
          else if (cli$present('Orbit')) Then
             pix_type = 4
             status = upm_get_float('Orbit',phi)
             write(asc_f,fmt='(f6.1)') phi
             c_line = c_line(:c_len) // '/TYPEPIX=ORBIT=' // asc_f
             c_len  = c_len + 21
          else
             pix_type = 1
             Big_Pixel = 6
             c_line = c_line(:c_len) // '/TYPEPIX=SKY=6'
             c_len  = c_len + 14 
          endif
          if (.not. status) call lib$signal(upm_invnum)
       end if
c*
c* Get galactic plane exclusion
c*
      if (status) then
          if (cli$present('GALAXY')) Then
              status = upm_get_float('GALAXY',galaxy)
              if (.not. status) call lib$signal(upm_invnum)
          else
              Galaxy = 0.
          end if
          write(asc_f,fmt='(f5.1)') galaxy
          c_line = c_line(:c_len) // '/GALAXY=' // asc_f
          c_len  = c_len + 13
      endif
c*
c* Get instrument quality thresholds (default=3, many yellow, no red flags)
c*
      if (status) then
          status = upm_get_longword('QUALITY.INSTRUMENT', Instr_Qual)
          status = upm_get_longword('QUALITY.ATTITUDE',   Attit_Qual)
          if (.not. status) call lib$signal(upm_invnum)
          write(ASC_I1,fmt='(i4)') instr_qual
          write(ASC_I2,fmt='(i4)') attit_qual
          c_line = c_line(:c_len) //
     .             '/QUAL=(IN=' // ASC_I1(4:4) // ',AT=' // ASC_I2(4:4) // ')'
          c_len  = c_len + 17
      endif
c*
c* Get orbital period (only if orbital pixels are used)
c*
      if ((pix_type .eq. 4) .and. status) then
          status = fut_orbital_period(jstart, R4_Orbital_PD)
          if (.not. status) call lib$signal(%val(status))

C*        Correction:  R4_period is in MINUTES 
C*                     R8_period is in 10^-7 SECONDS

          R8_Orbital_PD = 60.*R4_Orbital_PD * 1.E7	! 100nsec
          call AUT_DFloat2ADT(R8_Orbital_PD,period)
          call LIB$Subx(Time_Zero,period,period)
      endif
c*
c* Check for First Run, closeold, write, and fakeit qualifiers
c*
      if (status) then
          First_Run = cli$present('FIRST_RUN')
          if (first_run) then
              c_line = c_line(:c_len) // '/FIRST_RUN'
              c_len  = c_len + 10
          else
              c_line = c_line(:c_len) // '/NOFIRST_RUN'
              c_len  = c_len + 12
          endif
          Closeold  = cli$present('CLOSEOLD')
          if (closeold) then
              c_line = c_line(:c_len) // '/CLOSEOLD'
              c_len  = c_len + 9
          else
              c_line = c_line(:c_len) // '/NOCLOSEOLD'
              c_len  = c_len + 11
          endif
          Write = cli$present('WRITE')
          if (write) then
              c_line = c_line(:c_len) // '/WRITE'
              c_len  = c_len + 6
          else
              c_line = c_line(:c_len) // '/NOWRITE'
              c_len  = c_len + 8
          endif
          Fakeit = cli$present('FAKEIT')
          if (fakeit) then
              c_line = c_line(:c_len) // '/FAKEIT'
              c_len  = c_len + 7
              fakeit = 1
          else
              c_line = c_line(:c_len) // '/NOFAKEIT'
              c_len  = c_len + 9
              fakeit = 0
          endif
      endif
c*
c* Get the channels to process   
c*
      if (status) then
          status = cli$get_value ('CHANNEL', channel, length)
          if (.not. status) then
              call lib$signal(%val(status))
          else
              if (channel .eq. 'ALL') then
                  start_chan = 1
                  stop_chan  = 4
                  chan_int   = 1
              else if (channel .eq. 'RH') then
                  start_chan = 1
                  stop_chan  = 1
                  chan_int   = 1
              else if (channel .eq. 'RL') then
                  start_chan = 2
                  stop_chan  = 2
                  chan_int   = 1
              else if (channel .eq. 'LH') then
                  start_chan = 3
                  stop_chan  = 3
                  chan_int   = 1
              else if (channel .eq. 'LL') then
                  start_chan = 4
                  stop_chan  = 4
                  chan_int   = 1
              else if (channel .eq. 'RIGHT') then
                  start_chan = 1
                  stop_chan  = 2
                  chan_int   = 1
              else if (channel .eq. 'LEFT') then
                  start_chan = 3
                  stop_chan  = 4
                  chan_int   = 1
              else if (channel .eq. 'HIGH') then
                  start_chan = 1
                  stop_chan  = 3
                  chan_int   = 2
              else if (channel .eq. 'LOW') then
                  start_chan = 2
                  stop_chan  = 4
                  chan_int   = 2
              endif    
          endif
          c_line = c_line(:c_len) // '/CHAN=' // channel(:length)
          c_len  = c_len + 6 + length
       endif
c*
c* Get MTM mode to process
c*
      if (status) then
          num_mtm = 1
          c_line = c_line(:c_len) // '/MTM=('
          c_len  = c_len + 6
          status = cli$get_value('MTM', mtm(num_mtm), length)
          if ( .not. status ) call lib$signal(%val(status))
          if (mtm(num_mtm)(:3) .eq. 'ALL') then 
              c_line = c_line(:c_len) // 'ALL)' 
              c_len  = c_len + 4
          elseif (mtm(num_mtm)(:3) .ne. 'ALL') then
              c_line = c_line(:c_len) // mtm(num_mtm)(:length) 
              c_len  = c_len + length 
              do while (status .and. (num_mtm .le. 4))
                 num_mtm = num_mtm + 1
                 status = cli$get_value('MTM', mtm(num_mtm), length)
                 if ((.not. status) .and. (status .ne. %loc(cli$_absent))) then 
                      call lib$signal(%val(status))
                 elseif (status) then
                      c_line = c_line(:c_len) // ' ' // mtm(num_mtm)(:length) 
                      c_len  = c_len + 1 + length
                 endif
               enddo
               num_mtm = num_mtm -1
               if (status .eq. %loc(cli$_absent)) status = %loc(fss_normal)
               c_line = c_line(:c_len) // ')'
               c_len  = c_len + 1
           endif

           if (mtm(1)(:3) .eq. 'ALL') then
               mtmlist(1) = 1
               mtmlist(2) = 1
               mtmlist(3) = 1
               mtmlist(4) = 1
           else
               do i=1,num_mtm
                  if (mtm(i)(:2) .eq. 'SS') mtmlist(1) = 1
                  if (mtm(i)(:2) .eq. 'SF') mtmlist(2) = 1
                  if (mtm(i)(:2) .eq. 'LS') mtmlist(3) = 1
                  if (mtm(i)(:2) .eq. 'LF') mtmlist(4) = 1
               enddo
           endif

      endif
c*
c* Get pixel list to process
c*
      if (status .and. cli$present('PIXEL') ) then
          i = 1
          c_line = c_line(:c_len) // '/PIXEL = ( '
          c_len  = c_len + 11
          do while (status)
             status = upm_get_longword('PIXEL', plist(i))
             if ((.not. status) .and. (status .ne. absent_upm) ) then
                  call lib$signal(upm_invnum)
             elseif (status .ne. absent_upm) then
                  write(asc_i1,fmt='(i4)') plist(i)
                  c_line = c_line(:c_len) // asc_i1 // ' '
                  c_len  = c_len + 5
                  i = i + 1
             endif
           enddo
           c_line = c_line(:c_len) // ')'
           c_len  = c_len + 1
           if (status .eq. absent_upm) status = %loc(fss_normal)
           plist(i) = -9999

      else if (status .and. .not. cli$present('PIXEL') .and. 
     .		cli$present('PIXLIST') ) then
c*   Initilize the pixel list
	   do ilis=1,100 
		plists(ilis) = -1
	   enddo
c*   Read in the pixel list through as a namelist
	   
	   statusp = lib$get_lun( pixunit )
	   if ( .not. status ) call lib$signal( %val( status ) )

	   if ( statusp ) then
	     open ( pixunit, file='csdr$firas_out:pixlist.lis', status='OLD', 
     .		iostat=iosp )
	     if ( iosp .ne. 0 ) then
		call errsns( ios,,,,statusp)
		call lib$signal( %val( statusp ) )
		call lib$stop( %val( statusp ) )
	     endif
	     read ( unit=pixunit, NML=PIXLIST )
	     close ( pixunit )
	     ilis = 1
             c_line = c_line(:c_len) // '/PIXEL = ( '
             c_len  = c_len + 11
	     do while ( ( ilis .le. 99 ) .and. plists(ilis) .gt. -1 )
		i = ilis
		plist(i) = plists(i)
                write(asc_i1,fmt='(i4)') plist(i)
                c_line = c_line(:c_len) // asc_i1 // ' '
                c_len  = c_len + 5
		ilis = ilis + 1
	     enddo
c*   Take into account the case plists(100) > -1
	     if ( ( ilis .eq. 100 ) .and. plists(ilis) .gt. -1 ) then
		i = ilis
		plist(i) = plists(i)
                write(asc_i1,fmt='(i4)') plist(i)
                c_line = c_line(:c_len) // asc_i1 // ' '
                c_len  = c_len + 5
		ilis = ilis + 1
	     endif		
	   endif
c
c
c
	   if ( i .lt. 100 ) plist(i+1) = -9999
           c_line = c_line(:c_len) // ')'
           c_len  = c_len + 1

      else if (status .and. .not. cli$present('PIXEL') .and. 
     .		.not. cli$present('PIXLIST') ) then
           plist(1) = 9999
           c_line = c_line(:c_len) // '/PIXEL=(ALL)'
           c_len  = c_len + 12
      endif
c*
c* Get microprocessor mode list to process
c*
      if (status .and. cli$present('MPMODES') ) then
          i = 1
          c_line = c_line(:c_len) // '/MPMODES=('
          c_len  = c_len + 10
          do while (status)
             status = upm_get_longword('MPMODES', mpmodes(i))
             if ((.not. status) .and. (status .ne. absent_upm) ) then
                  call lib$signal(upm_invnum)
             elseif (status .ne. absent_upm) then
                  write(asc_i1,fmt='(i4)') mpmodes(i)
                  c_line = c_line(:c_len) // asc_i1 // ' '
                  c_len  = c_len + 5
                 i = i + 1
             endif
           enddo
           c_line = c_line(:c_len) // ')'
           c_len  = c_len + 1
           if (status .eq. absent_upm) status = %loc(fss_normal)
           mpmodes(i) = -9999
      else if (status .and. .not. cli$present('MPMODES')) then
           mpmodes(1) = 9999
           c_line = c_line(:c_len) // '/MPMODES=(ALL)'
           c_len  = c_len + 14
      endif
c*
c* Get time for oldest data to process
c*
      if (status) then
          status = cli$get_value('OLDEST', oldest, LENGTH)
          if (.not. status) call lib$signal(%val(status))
          c_line = c_line(:c_len) // '/OLDEST=' // oldest(:length)
          c_len  = c_len + 8 + length
      endif      
c
c Was input skymap defaulted ?
c
      if (status) then
          if (cli$present('INSKYMAP')) then
              in_def = .false.
          else
              in_def = .true.
         endif
      endif
c
c Get input skymap file name
c
      if (status) then
          if (cli$present('INSKYMAP')) then
              status = cli$get_value('INSKYMAP', inskymap, length)
              if (.not. status ) then
                  call lib$signal(%val(status))
              else
                  c_line = c_line(:c_len) // '/INSKYMAP=' // inskymap(:length)
                  c_len  = c_len + 10 + length
              endif
          else
              inskymap = ' '
          endif
      endif          
c
c Was output skymap defaulted ?
c
      if (status) then
          if (cli$present('OUTSKYMAP')) then
              out_def = .false.
          else
              out_def = .true.
         endif
      endif
c
c Get output skymap file name
c
      if (status .and. write) then
          if (cli$present('OUTSKYMAP')) then
              status = cli$get_value('OUTSKYMAP', outskymap, length)
              if (.not. status) then
                  call lib$signal(%val(status))
              else
                  c_line = c_line(:c_len) // '/OUTSKYMAP=' // outskymap(:length)
                  c_len  = c_len + 11 + length
              endif
          else
              outskymap = ' '
          endif
      elseif (status .and. .not. write) then
          outskymap = ' '
      endif          
c
c Was eng file defaulted ?
c
      if (status) then
          if (cli$present('ENGFILE')) then 
              eng_def = .false.
          else
              eng_def = .true.
         endif
      endif
c
c Get engineering file name
c
      if (status) then
          if (cli$present('ENGFILE')) then
              status = cli$get_value('ENGFILE', engfile, length)
              if ( .not. status ) call lib$signal(%val(status))
              c_line = c_line(:c_len) // '/ENGFILE=' // engfile(:length)
              c_len  = c_len + 9 + length
          else
              engfile = ' '
          endif
      endif          
c
c Get RSE file name. Since there is no default RSE, a file name 
c must be supplied. The CLD requires a file name with /RSE.
c
      if (status) then
          c_line = c_line(:c_len) // '/RSE=' 
          c_len  = c_len + 5 

          rse_stat = cli$get_value('RSE', rsefile, length)

          if (rse_stat .and. (length .gt. 0)) then
              c_line = c_line(:c_len) // rsefile(:length)
              c_len  = c_len + length
          else if (.not. rse_stat .or. (length .eq. 0)) then
              rsefile = ' '
              c_line = c_line(:c_len) // 'NORSE'
              c_len  = c_len + 5
          endif
      endif          
c*
c* Get ICAL temp tolerance 
c*
      if (status) then
          status = upm_get_float('TEMPTOL',temptol)
          if (.not. status) then
               call lib$signal(upm_invnum)
          else
               write(atemptol,fmt='(f5.1)') temptol
               c_line = c_line(:c_len) // '/TEMPTOL=' // atemptol
               c_len  = c_len + 14
               temptol = temptol / 1000.0  ! Units are millikelvin
          endif
      endif
c*
c* Get minimum Sun angle allowed
c*
      if (status) then
          status = upm_get_float('MINSUN',minsun)
          if (.not. status) then
               call lib$signal(upm_invnum)
          else
               write(aminsun,fmt='(f6.2)') minsun
               c_line = c_line(:c_len) // '/MINSUN=' // aminsun
               c_len  = c_len + 14
          endif
      endif
c*
c* Get maximum Sun angle allowed
c*
      if (status) then
          status = upm_get_float('MAXSUN',maxsun)
          if (.not. status) then
               call lib$signal(upm_invnum)
          else
               write(amaxsun,fmt='(f6.2)') maxsun
               c_line = c_line(:c_len) // '/MAXSUN=' // amaxsun
               c_len  = c_len + 14
          endif
      endif
c*
c* Get minimum Moon angle allowed
c*
      if (status) then
          status = upm_get_float('MINMOON',minmoon)
          if (.not. status) then
               call lib$signal(upm_invnum)
          else
               write(aminmoon,fmt='(f6.2)') minmoon
               c_line = c_line(:c_len) // '/MINMOON=' // aminmoon
               c_len  = c_len + 15
          endif
      endif
c*
c* Get maximum Moon angle allowed
c*
      if (status) then
          status = upm_get_float('MAXMOON',maxmoon)
          if (.not. status) then
               call lib$signal(upm_invnum)
          else
               write(amaxmoon,fmt='(f6.2)') maxmoon
               c_line = c_line(:c_len) // '/MAXMOON=' // amaxmoon
               c_len  = c_len + 15
          endif
      endif
c*
c* Get minimum Earth limb angle allowed
c*
      if (status) then
          status = upm_get_float('MINEARTH',minearth)
          if (.not. status) then
               call lib$signal(upm_invnum)
          else
               write(aminearth,fmt='(f6.2)') minearth
               c_line = c_line(:c_len) // '/MINEARTH=' // aminearth
               c_len  = c_len + 16
          endif
      endif
c*
c* Get maximum Earth limb angle allowed
c*
      if (status) then
          status = upm_get_float('MAXEARTH',maxearth)
          if (.not. status) then
               call lib$signal(upm_invnum)
          else
               write(amaxearth,fmt='(f6.2)') maxearth
               c_line = c_line(:c_len) // '/MAXEARTH=' // amaxearth
               c_len  = c_len + 16
          endif
      endif
c*
c* Get ICAL temperature bins. 
c*
      bins_defaulted = cli$present('BINS') .eq. %loc(cli$_defaulted) 

      if (status) then
         if (bins_defaulted) then
             num_bins = 3
             bins(1)  = 2.747
             bins(2)  = 2.758
             bins(3)  = 2.771
             c_line = c_line(:c_len) // '/BINS=(2.747 2.758 2.771)'
             c_len  = c_len + 25
         else
             i        = 1
             num_bins = 0
             c_line = c_line(:c_len) // '/BINS=('
             c_len  = c_len + 7
             do while (status)
                status = upm_get_float('BINS', bins(i))
                if ((status .eq. absent_upm) .and. (i .eq. 1)) then
                    num_bins = 3
                    bins(1)  = 2.747
                    bins(2)  = 2.758
                    bins(3)  = 2.771
                    c_line = c_line(:c_len) // '2.747 2.758 2.771'
                    c_len  = c_len + 17
                elseif ((.not. status) .and. (status .ne. absent_upm)) then
                    call lib$signal(upm_invnum)
                elseif (status .ne. absent_upm) then
                    write(asc_f,fmt='(f6.3)') bins(i)
                    c_line = c_line(:c_len) // asc_f
                    c_len  = c_len + 6
                    i = i + 1
                    num_bins = num_bins + 1
                endif
             enddo
             c_line = c_line(:c_len) // ')'
             c_len  = c_len + 1
             if (status .eq. absent_upm) status = %loc(fss_normal)
         endif
      endif
c
c
c Check temp bins for overlap. If there is overlap, signal error.
c
      overlap = .false.
      do i = 1,num_bins-1
         rtemp = bins(i) + (2.0 * temptol)
         if ( rtemp .gt. bins(i+1) ) overlap = .true.
      enddo

      if (overlap) then
          call lib$signal(fss_bin_overlap)              
          status = %loc(fss_bin_overlap)
      endif   
c
c* Get Dihedral temperature bins. 
c*
      bins_defaulted = cli$present('DBINS') .eq. %loc(cli$_defaulted) 

      if (status) then
         if (bins_defaulted) then
             num_dbins = 3
             dbins(1)  = 3.5
             dbins(2)  = 4.5
             dbins(3)  = 5.5
             c_line = c_line(:c_len) // '/DBINS=(3.50 4.50 5.50'
             c_len  = c_len + 22
         else
             i        = 1
             num_dbins = 0
             c_line = c_line(:c_len) // '/DBINS=('
             c_len  = c_len + 8
             do while (status)
                status = upm_get_float('DBINS', dbins(i))
                if ((status .eq. absent_upm) .and. (i .eq. 1)) then
                    num_dbins = 3
                    dbins(1)  = 3.5
                    dbins(2)  = 4.5
                    dbins(3)  = 5.5
                    c_line = c_line(:c_len) // '3.50 4.50 5.50'
                    c_len  = c_len + 14
                elseif ((.not. status) .and. (status .ne. absent_upm)) then
                    call lib$signal(upm_invnum)
                elseif (status .ne. absent_upm) then
                    if (i .gt. 1) then
                       write(asc_f,fmt='(f5.2)') dbins(i)
                       c_line = c_line(:c_len) // asc_f
                       c_len  = c_len + 5
                       i = i + 1
                       num_dbins = num_dbins + 1
                    endif
                    if (i .eq. 1) then
                       write(asc_f,fmt='(f4.2)') dbins(i)
                       c_line = c_line(:c_len) // asc_f
                       c_len  = c_len + 4
                       i = i + 1
                       num_dbins = num_dbins + 1
                    endif
                endif
             enddo
             c_line = c_line(:c_len) // ')'
             c_len  = c_len + 1
             if (status .eq. absent_upm) status = %loc(fss_normal)
         endif
      endif
c
c
c Check that dihedral temp bins are sorted. If not, signal error.
c
      dihed_sort = .false.
      do i = 1,num_dbins-1
         if ( dbins(i) .ge. dbins(i+1) ) dihed_sort = .true.
      enddo

      if (dihed_sort) then
          call lib$signal(fss_dihed_sort)              
          status = %loc(fss_dihed_sort)
      endif   
c
c Break up command line into one-line segments.
c First initialize command array c_array to blanks.
c
      call str$trim(c_line, c_line, c_len)
      
      do i = 1,30
         c_array(i) = ' '
      enddo

c Find location of all slashes in command line.        

      num_slash = 0
      do i = 1,c_len
         if (c_line(i:i) .eq. '/') then
             num_slash = num_slash + 1
             slash(num_slash) = i
         endif
      enddo

c Put in a slash at the end to make processing easier    

      num_slash = num_slash + 1
      slash(num_slash) = c_len + 1
c
c Go through command line, putting a substring in array if 
c it is >= 75 chars long and the next slash would make the 
c line too long. Break a qualifier if it is longer than 75 chars.
c
      line_start = 1
      cnum       = 0
      qual_len   = 0
      i = 2
      do while (i .le. num_slash)
         if ((slash(i) .lt. (line_start+line_len)).and.(i .lt. num_slash)) then
             continue
         elseif (i .eq. num_slash) then
             if (slash(i) .le. line_start+line_len) then
                 cnum = cnum + 1
                 c_array(cnum) = c_line(line_start:slash(i)-1)
             else
                 if ( c_line(slash(i-1):slash(i-1)+4) .eq. '/BINS' ) then
                    cnum = cnum + 1
                    c_array(cnum) = c_line(line_start:slash(i-1)-1)
                    line_start = slash(i-1)
                    use_line_len = pix_line_len
                 endif
                 do while (slash(i) .gt. (line_start+1) )
                    j = min(slash(i)-1,line_start+line_len)
                    cnum = cnum + 1
                    c_array(cnum) = c_line(line_start:j)
                    line_start = j + 1
                 enddo
             endif
         elseif ((slash(i).le.(slash(i-1)+line_len)).and.(i .lt. num_slash))then
             cnum = cnum + 1
             c_array(cnum) = c_line(line_start:slash(i-1)-1)
             line_start = slash(i-1) 
         elseif ((slash(i).gt.(slash(i-1)+line_len)).and.(i.lt.num_slash))then
             if ( c_line(slash(i-1):slash(i-1)+5) .eq. '/PIXEL' ) then
                cnum = cnum + 1
                c_array(cnum) = c_line(line_start:slash(i-1)-1)
                line_start = slash(i-1)
                use_line_len = pix_line_len
             else
                use_line_len = line_len
             endif                    ! End break in the com line for pixel
             do while (slash(i) .gt. (line_start+1) )
                j = min(slash(i)-1,line_start+use_line_len)
                cnum = cnum + 1
                c_array(cnum) = c_line(line_start:j)
                line_start = j + 1
             enddo
         endif
         i = i + 1
      enddo             

      fss_get_quals = status

      return
      end
