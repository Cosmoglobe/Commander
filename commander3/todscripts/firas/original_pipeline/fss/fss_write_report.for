      integer*4 function fss_write_report(num_groups, bad_groups, good_ifgs,
     .                                    bad_ifgs, pix_array, mode_array, 
     .                                    bin_array, dbin_array, size_array,  
     .                                    minimum, temptol, detail, chan, 
     .                                    lu_rep,in_ext,out_ext,eng_ext,frun,
     .                                    write)

C------------------------------------------------------------------------
C    PURPOSE: This function writes the summary information for one
C             channel of the FSS run to a report file. 
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION:  status=fss_write_report(num_groups, bad_groups, good_ifgs,
C     .                                   bad_ifgs, pix_array, mode_array,
C     .                                   bin_array, dbin_array, size_array,
C     .                                   minimum, temptol, detail, chan, 
C     .                                   lu_rep,in_ext,out_ext,eng_ext,frun,
C     .                                   write)
C
C
C    INPUT PARAMETERS:   I*4   num_groups              total coadd grps 
C                        I*4   bad_groups              coadd grps < minimum
C                        I*4   good_ifgs               number coaddable ifgs
C                        I*4   bad_ifgs                number uncoaddable ifgs
C                        I*4   pix_array (max_groups)  pixel numbers for coadds
C                        I*4   mode_array(max_groups)  modes for coadds
C                        I*4   size_array(max_groups)  size of coadd groups
C                        I*4   bin_array (max_groups)  ICAL temp bin of coadds
C                        I*4   dbin_array (max_groups) dihedral temp bin of coadds
C                        I*4   minimum(4)              minimum size of coadd
C                        R*4   temptol                 Temp tolerance
C                        c*10  detail                  Level of detail 
C                        I*4   chan                    current channel to report
C                        I*4   lu_rep                  logical unit for report
c                        C*32  in_ext                  input skymap extension
c                        C*32  out_ext                 output skymap extension
c                        C*32  eng_ext                 eng file extension
c                        I*4   frun                    first run
c                        I*4   write                   write/nowrite flag
c
C
C    OUTPUT PARAMETERS:  None. 
C
C
C    SUBROUTINES CALLED:    cut_register_version
C                           cut_display_banner
C                           cut_translate_archive_id
C                           str$upcase
C
C    COMMON VARIABLES USED: None.
C 
C    INCLUDE FILES:         $jpidef
C                           fss_include
C
C----------------------------------------------------------------------
C
C Changes:
C
C       SPR 12197 - Modifications to implement /DBINS qualifier in place
C                   of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C----------------------------------------------------------------------
C                  PDL FOR fss_write_report
C
C                  D. BOULER, JAN, 1991
C
C fss_write_report = %loc(fss_normal)         
C status           = %loc(fss_normal)
C
C Convert binary time tolerance to ascii 
C
c if (first channel) then
c     Write input and output skymap and engineering names to report
c endif
c
C write summary report lines for channel
C increment lines
C
C if (status .and. (detail level is "DETAIL")) then
C
C     grps = 0
C     do while (grps .le. num_grps(chan))
C
C        if ((lines .ge. max_lines) .and. (grps .lt. num_grps(chan)) then
C             write page header with channel number
C             write header for detail line
C             increment lines
C        elseif (grps .eq. 0) then
C             write header for detail line
C             increment lines
C        endif
C
C        write detail line
C        increment grps with number of groups per detail line
C        increment lines with number of lines per detail line
C
C        if ((grps .ge. num_groups) .and. (channel not last)) then
C            write line at top of page for next channel
C        endif
C
C       enddo    !  while (grps .le. num_grps(chan))
C  endif         !  (status .and. (detail level is "DETAIL"))
c
C fss_write_report = status 
C
C return
C end (pdl)
C*****************************************************************************
C
      implicit none
c
c Include files:
c
      include  '($jpidef)'
      include  '(fss_include)'
c
c Declare externals
c
      external fss_normal
      external cut_translate_archive_id
c
c Function declarations:
c
      integer*4 cut_translate_archive_id
      integer*4 cut_register_version
      integer*4 cut_display_banner
      integer*4 lib$get_lun
C
C Variable declarations
C
      integer*4 i,j,k,l                    ! Loop counters
      integer*4 chan                       ! Current channel
      integer*4 status                     ! General status return
      integer*4 ios                        ! i/o status return
      integer*4 good_ifgs                  ! Number of good ifgs/channel
      integer*4 bad_ifgs                   ! Number of bad ifgs/channel
      integer*4 num_groups                 ! Total number coadd groups/channel
      integer*4 bad_groups                 ! Number small coadd groups/channel
      integer*4 pix_array (fss_max_groups) ! Pixel numbers for coadd groups
      integer*4 size_array(fss_max_groups) ! Size of coadd groups
      integer*4 mode_array(fss_max_groups) ! Modes for coadd groups
      integer*4 bin_array (fss_max_groups) ! ICAL Temp bins for coadd groups
      integer*4 dbin_array (fss_max_groups) ! Dihedral Temp bins for coadd groups
      integer*4 lu_rep                     ! Logical unit for report
      integer*4 len                        ! Length of report file name
      integer*4 temp_len                   ! Length of temporary logical name 
      integer*4 trans_len(8)               ! Length of translated logical name
      integer*4 today(2)                   ! Today's date (adt format)
      integer*4 lines                      ! Number of lines/page written
      integer*4 max_lines      /55/        ! Max number of lines/page
      integer*4 grps_per_line  /6/         ! Groups per detail line
      integer*4 grps_left                  ! Number of remaining groups
      integer*4 grps_to_print              ! Number of groups to print on line
      integer*4 minimum(4)                 ! Minimum number of ifgs/group
      integer*4 grps                       ! Current group number
      integer*4 frun                       ! FIRST_RUN true?
      integer*4 write                      ! WRITE specified ?

      real*4    temptol                ! Temp tolerance
      
      character*2    achan(4)          ! Ascii version of channel (rh, etc.)
      data           achan                    /'RH','RL','LH','LL'/

      character*10   detail            ! Level of detail for report

      character*14   gmt               ! Today's data (gmt format)

      character*72   trans_logs(8)     ! Translated logicals 
      character*72   temp_log          ! Temporary logical name
      character*72   in_logs(8)        ! Input logicals (CSDR$FIRAS_IN, eg)

      character*128  inskymap          ! Input skymap name
      character*128  outskymap         ! Output skymap name
      character*128  engfile           ! Engineering file name
      character*128  in_ext            ! Input  skymap extension
      character*128  out_ext           ! Output skymap extension
      character*128  eng_ext           ! Eng file extension

      logical*4      first  /.true./   ! First pass through?          

      save           first             ! Save value of first
c*
c******  Begin code *******************************************
c*
      fss_write_report = %loc(fss_normal)         
      status           = %loc(fss_normal)         
c
c If first channel, write input and output skymap names to report.
c
       if (first) then

           call str$trim(in_ext, in_ext)
           if (frun) then
               inskymap  = 'NONE          : FIRST_RUN specified'
           else
               inskymap  = 'CSDR$FIRAS_IN :FSS_SSSKY_XX.' //  in_ext
           endif

           call str$trim(out_ext, out_ext)
           if (write) then
               outskymap = 'CSDR$FIRAS_OUT:FSS_SSSKY_XX.' // out_ext
           else
               outskymap = 'NONE          : NOWRITE   specified'
           endif
          
           call str$trim(eng_ext, eng_ext)
           engfile =   'CSDR$FIRAS_RAW:FDQ_ENG.' // eng_ext

           write(lu_rep,*) ' '
           write(lu_rep,*) 'Input  skymap: ',inskymap(:60)
           write(lu_rep,*) 'Output skymap: ',outskymap(:60)
           write(lu_rep,*) 'Eng file     : ',engfile(:60)
           write(lu_rep,*) ' '
       endif 
c
c Write summary lines for channel
c
      if (first) then
          lines = 13
      else
          lines = 9
      endif

      write(lu_rep,20)  achan(chan)
      write(lu_rep,40)  minimum(chan), temptol, 
     .                 (good_ifgs + bad_ifgs), good_ifgs, bad_ifgs, 
     .                  num_groups, (num_groups - bad_groups), bad_groups

20    format(1x,'****************** Channel ', a, ' ******************',/)
30    format('1','****************** Channel ', a, ' ******************',/)
40    format(1x, ' Minimum = ',I4, '  Temp Tolerance = ', f7.4, //,
     .       1x, I9, ' total interferograms input for processing', /,
     .       1x, I9, ' interferograms went into groups >= minimum', /,
     .       1x, I9, ' interferograms went into groups <  minimum', /,
     .       1x, I9, ' total coadd groups were formed', /,
     .       1x, I9, ' coadd groups with >= minimum were formed', /,
     .       1x, I9, ' coadd groups with <  minimum were formed', /)

      lines = lines + 9     
c
c Write detail lines if detail format is requested 
c
      if (detail .eq. 'FULL') then 

          grps = 0

          do while (grps .lt. num_groups)

             if ((lines .ge. max_lines).and.(grps .le. num_groups)) then
                  write(lu_rep,30) achan(chan)
                  write(lu_rep,50) 
50                format(x,6('PIXL S  B NUM '))
                  lines = 3
             elseif (grps .eq. 0) then
                  write(lu_rep,50) 
                  lines = lines + 1
             endif    

             grps_left     = num_groups - grps
             grps_to_print = grps + min(grps_per_line,grps_left) 

             write(lu_rep,60) (pix_array(l),mode_array(l),bin_array(l),dbin_array(l),
     .                         size_array(l), l = grps+1,grps_to_print)
60           format(x,<grps_per_line>(i4,x,i1,x,i2,x,i2,x,i3,x))

             grps  = grps  + grps_per_line
             lines = lines + 1

           enddo         ! while (grps .le. num_groups
      endif              ! (detail .eq. "FULL") 

      if (first) first = .false.

      fss_write_report = status 

      return
      end
