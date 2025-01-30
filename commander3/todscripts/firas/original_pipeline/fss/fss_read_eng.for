      integer*4 function fss_read_eng( jstart, jstop, rse, num_eng, xcal,
     .                                 ical, refh, skyh, dihd, bolo, sci_times,
     .                                 grtwt, grttrans,start_chan, stop_chan,
     .                                 chan_int, lu_rep, instr_qual,
     .                                 attit_qual, mpmodes, mtmlist, oldest,
     .                                 bins, num_bins, temptol, max_dihed,
     .                                 eng_ext, fakeit, report)

C------------------------------------------------------------------------
C    PURPOSE: Reads the engineering data. Applies an rse to the data and 
C             returns only the data which satisfy the rse and command
C             line parameters. If no RSE is specified, reads all eng data 
C             in the time range.
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C
C    INVOCATION:  status = fss_read_eng(jstart, jstop, rse, num_eng, xcal, ical,
C                                       refh, skyh, dihd,bolo,sci_times,grtwt,
C                                       grttrans,start_chan,stop_chan,chan_int,
C                                       lu_rep,instr_qual, attit_qual, mpmodes,
C                                       mtmlist,oldest,bins,num_bins, temptol, 
C                                       max_dihed, eng_ext, fakeit, report)
C
C    INPUT PARAMETERS:    c*14  jstart                 Start time of processing
C                         c*14  jstop                  End   time of processing
C                         c*128 rse(16)                Record selection exp
C                         rec   grtwt                  GRT weights
C                         rec   grttrans               GRT transition temps
c                         i*4   start_chan             start channel
c                         i*4   stop_chan              stop  channel
c                         i*4   chan_int               channel interval
C                         i*4   lu_rep                 report file LU
C                         I*4   instr_qual             instrument quality limit
C                         i*4   attit_qual             attitude quality limit
C                         i*4   mpmodes(*)             micro modes to process
C                         i*4   mtmlist(4)             MTM modes to process
C                         c*14  oldest                 Earliest acceptable time
C                         R*4   bins(*)                ICAL temp bins
C                         i*4   num_bins               Number of ICAL temp bins
C                         r*4   temptol                temp tolerance
C                         r*4   max_dihed              Maximum dihedral temp
C                         c*128 eng_ext                Input eng file extension
C                         i*4   fakeit                 Is fakeit data allowed?
C                         I*4   report                 Report requested?
C
C    OUTPUT PARAMETERS:   i*4   num_eng(4)              Number of eng recs/chan
C                         r*4   xcal(4,max_buff)        xcal temps
C                         r*4   ical(4,max_buff)        ical temps
C                         r*4   refh(4,max_buff)        ref horn temps
C                         r*4   skyh(4,max_buff)        sky horn temps
C                         r*4   dihd(4,max_buff)        dihedral temps
C                         r*4   bolo(4,max_buff)        bolometer temps
C                         i*4   sci_times(2,4,max_buff) Science binary times
C
C    SUBROUTINES CALLED:  CT_QUERY_ARCV
C                         CT_QUERY_GET
C                         CT_OPEN_QUERY_NONRAW
C                         FUT_TEMP_LIST
C
C    COMMON VARIABLES USED: None.
C 
C    INCLUDE FILES:   ct$library:ctuser.inc
C                     fut_params
C                     fss_include
C                     cct_query_catalog_record
C
c---------------------------------------------------------------------
c                   PDL for FSS_READ_ENG
c                   D. Bouler, Oct, 1990
c
C  Get engineering LU 
C
C  Find whether RSE was specified on command line
C
C  If (eng file name not specified) then
C      build engineering file name 
C  else
C      use eng file specified on command line
C      get JSTART and JSTART for file from catalog; override command line
C  endif
C
C  if (no RSE specified) then
C      open engineering archive for sequential access
C  else
c      open eng file through query path
C      Call CT_QUERY_ARCV to apply RSE to input data
C  endif
C
C  Do while ((not end of data) and. (number of recs .lt. maximum))
C     
C     Read engineering record
C     if (COBETRIEVE status) then
C         do i = start_chan, stop_chan, chan_interval
C            if ((sci_time(i) .le. start) .and. (sci_time(i) .ge. stop)) then
C                 num_eng(i) = num_eng(i) + 1
C                 if (num_eng(i) .eq. maximum) signal max data read in 
C
C                 sci_times(1,i,num_eng(i)) = sci_time(1,i)
C                 sci_times(2,i,num_eng(i)) = sci_time(2,i)
C
C                 call fut_temp_list to get xcal, ical, ref horn, 
C                 sky horn, and dihedral  temps
C
C                 Store xcal,ical,ref horn,sky horn,dihedral temperatures
C
C            endif       !  {science time OK} 
C          enddo         !  {channels}
C       endif            !  {COBETRIEVE status OK}
C   enddo                !  {while not end of data and number of recs lt MAX}
C
C   Write summary of total records read and rejected to report file
C
C   if (no eng recs found) then
C       signal error
C       status = %loc(fss_noengrecs)
C   endif
C
C   call CT_CLOSE_ARCV for eng archive
C   if (COBETRIEVE status not OK) then
C       signal error
C       status = %loc(fss_badengclose)
C   endif
C
C   fss_read_eng = status
C
C   return
C   end (pdl)
C----------------------------------------------------------------------
C
C Changes:
C
C   D. Bouler, 12/19/91, SPR 9364 -- Increase number of eng records that
C   can be read.
C
C
C   SPR 12197 - Modifications to implement /DBINS qualifier in place
C               of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C----------------------------------------------------------------------

      implicit none
c
c Include files
c
      Include 'CT$LIBRARY:CTUser.Inc'
      Include '(FUT_Params)'
      Include '(fss_include)'
c
c FUNCTION DECLARATIONS
c
      integer*4      lib$get_lun
      integer*4      sys$asctim
      integer*4      fut_temp_list

      logical*2      time_lt
      logical*2      time_gt
      logical*2      time_le
      logical*2      time_ge
c*
c* Declare externals
c*
      external       ct_connect_query
      external       ct_connect_read
      external       ct_query_arcv
      external       ct_binary_to_gmt
      external       fss_badkeyread
      external       fss_normal
      external       fss_badengopen
      external       fss_badengclose
      external       fss_badqget
      external       fss_badseqrd
      external       fss_badquery
      external       fss_noengrecs
      external       fss_allengrej
      external       fss_nocatrecs
      external       fss_maxdata
C
C Variable declarations
C
      Integer*2	     ct_status(20)         ! COBETRIEVE status return

      integer*4      start(2)              ! ADT start time
      integer*4      stop(2)               ! ADT stop time
      integer*4      sci_times(2,4,fss_max_buff) ! Sci rec times
      integer*4      i,j,k,l,m             ! Counters
      integer*4      temp_stat             ! Return status for checking temp
      integer*4      status                ! General status return
      integer*4      lu_eng                ! LU for engineering
      integer*4      ios                   ! I/O status return
      integer*4      adt(2)                ! Today's date in binary format
      integer*4      eng_ctr               ! Counter for eng recs 
      integer*4      num_eng(4)            ! Number of eng records/channel 
      integer*4      start_chan            ! Start channel
      integer*4      stop_chan             ! Stop  channel
      integer*4      chan_int              ! Channel interval
      integer*4      curr_time(2)          ! Current time and date
      integer*4      combswt    /0/        ! Combine B side only
      integer*4      out_time(4)           ! Recs rejected because of bad time
      integer*4      out_notime(4)         ! Recs rejected because no sci time
      integer*4      out_qual(4)           ! Counter for bad quality records
      integer*4      out_xcal(4)           ! Number of recs with xcal in cal
      integer*4      out_mode(4)           ! number of recs not in mode list
      integer*4      out_mtm(4)            ! number of recs not in mtm list
      integer*4      out_bin(4)            ! number of recs not in temp bin
      integer*4      out_fake(4)           ! number of recs with bad fakeit
      integer*4      out_dihed(4)          ! number of recs with high dihedral
      integer*4      all_eng(4)            ! Total eng recs read
      integer*4      lu_rep                ! Report file LU
      integer*4      stop_lu               ! Max LU to print info to
      integer*4      lu_int                ! Interval between LUs
      integer*4      tot_read              ! Total eng recs read
      integer*4      tot_good              ! Total eng recs not rejected
      integer*4      lu                    ! Generic LU
      integer*4      instr_qual            ! Instrument quality limit
      integer*4      attit_qual            ! Attitude quality limit
      integer*4      mpmodes(*)            ! Microprocessor modes to process
      integer*4      mtmlist(4)            ! MTM modes to process
      integer*4      earliest(2)           ! Earliest acceptable time (binary)
      integer*4      xcal_pos1             ! XCAL position 
      integer*4      xcal_pos2             ! XCAL position 
      integer*4      num_modes             ! Number of micro modes specified
      integer*4      sci_mode              ! Micro mode for record
      integer*4      mtm                   ! MTM mode for record
      integer*4      iqual                 ! Instrument quality for record
      integer*4      aqual                 ! Attitude quality for record
      integer*4      speed                 ! MTM speed
      integer*4      length                ! MTM length
      integer*4      num_bins              ! Number of ICAL temp bins
      integer*4      eng_len               ! Length of eng file name
      integer*4      period                ! Position of "." in file name
      integer*4      fakeit                ! Was fakeit allowed?
      integer*4      eng_fakeit            ! Fakeit setting for current record
      integer*4      report                ! Was report requested?
   
      real*4         xcal(4,fss_max_buff)  ! xcal temps
      real*4         ical(4,fss_max_buff)  ! ical temps
      real*4         refh(4,fss_max_buff)  ! ref horn temps
      real*4         skyh(4,fss_max_buff)  ! sky horn temps
      real*4         dihd(4,fss_max_buff)  ! sky horn temps
      real*4         bolo(4,fss_max_buff)  ! bolometer temps
      real*4         outtemps(10)          ! temp array for fut_temp_list
      real*4         outsigs(10)           ! output temp sigmas
      real*4         bins(*)               ! ICAL temp bins
      real*4         temptol               ! ICAL temp tolerance
      real*4         ical_temp             ! ICAL temp for record
      real*4         max_dihed             ! Maximum dihedral temp

      logical*4      mode_in_list          ! Is mode in list?
      logical*4      check_modes           ! Are sci modes being checked ?
      logical*4      mtm_in_list           ! Is mtm in list?
      logical*4      check_mtm             ! Are mtm modes being checked ?
      logical*4      singlifg  /.true./    ! Switch for single-ifg 
      logical*4      in_bin                ! Is ical temp in bin ?
      logical*4      eng_spec              ! Was eng file specified?
      logical*4      in_time_range         ! Is record in time range?

      character*14   jstart                ! Start time (gmt)
      character*14   jstop                 ! Stop time  (gmt)
      character*14   oldest                ! Earliest acceptable time
      character*14   grt_name(10)          ! GRT names

      character*30   gmt                   ! GMT time      
      character*30   time_str              ! Current date and time

      character*128  eng_ext               ! Name of engineering file ext
      character*128  engfile               ! Name of engineering file
      character*128  rse(16)               ! Record selection expression

      character*2    chan_id(4)/'RH','RL','LH','LL'/     ! Channel identifiers

c
c Declare records
c                                             
      dictionary 'fdq_eng'
      dictionary 'fex_grtrawwt'
      dictionary 'fex_grttrans'
      dictionary 'fut_enganlg'

      record /fdq_eng/          eng_rec    ! Engineering record
      record /fex_grtrawwt/     grtwt      ! GRT weights
      record /fex_grttrans/     grttrans   ! GRT transition temps
      record /fut_enganlg/      rawsigs    ! Raw sigmas 
c
c***  Begin code *******************************************
c
      fss_read_eng = %loc(fss_normal)  

      do i = 1,4
         out_xcal(i)   = 0
         out_qual(i)   = 0
         out_mode(i)   = 0
         out_mtm(i)    = 0
         out_bin(i)    = 0
         out_time(i)   = 0
         out_notime(i) = 0
         out_fake(i)   = 0
         out_dihed(i)  = 0
         all_eng(i)    = 0       
      enddo

      grt_name(01) = 'XCAL'
      grt_name(02) = 'ICAL'
      grt_name(03) = 'SKY HORN'
      grt_name(04) = 'REF HORN'
      grt_name(05) = 'DIHEDRAL'
      grt_name(06) = 'STRUCTURE'
      grt_name(07) = 'RH BOLOMETER'
      grt_name(08) = 'RL BOLOMETER'
      grt_name(09) = 'LH BOLOMETER'
      grt_name(10) = 'LL BOLOMETER'

      if (eng_ext .eq. ' ') then
          eng_spec = .false.
      else
          eng_spec = .true.
      endif

      if (mpmodes(1) .ne. 9999) then
          num_modes = 0
          do while (mpmodes(num_modes+1) .ne. -9999)  
             num_modes = num_modes + 1
          enddo
          check_modes  = .true.
      else
          check_modes  = .false.
          mode_in_list = .true.
      endif

      if ((mtmlist(1) .eq. 1) .and. (mtmlist(2) .eq. 1) .and.
     .    (mtmlist(3) .eq. 1) .and. (mtmlist(4) .eq. 1)) then
           check_mtm   = .false.
           mtm_in_list = .true.
      else
           check_mtm   = .true.
      endif

      call ct_gmt_to_binary(oldest, earliest)
c
c Build eng file name
c
      call str$trim(eng_ext, eng_ext, eng_len)

      if (eng_ext .eq. ' ') then
          engfile = 'CSDR$FIRAS_RAW:FDQ_ENG/'//jstart//';'//jstop// ';'
      else
          period = index(eng_ext, '.')
          if (period .gt. 0) then
              eng_ext = eng_ext(period+1:eng_len)
              call str$trim(eng_ext, eng_ext, eng_len)
              engfile = 'CSDR$FIRAS_RAW:FDQ_ENG.'//eng_ext(:eng_len)
          else
              engfile = 'CSDR$FIRAS_RAW:FDQ_ENG.' // eng_ext(:eng_len)
          endif
       endif
c
c If eng file specified, get jstart and jstop from file name.
c
        if (eng_spec) then
            jstart = eng_ext(4:8) // '000000000'
            jstop  = eng_ext(4:8) // '235959999'
        endif
c
c OPEN Engineering Archive
c
      status = lib$get_lun (LU_eng)
      if (.not. status) call lib$signal(%val(status))

      if (status .and. (rse(1) .ne. ' ') ) then 
          open ( unit     = LU_eng,
     .           file     = engfile,
     .           iostat   = ios,
     .           status   = 'OLD',
     .           useropen = CT_CONNECT_QUERY)
 
          if (ios .ne. 0) then
              status = %loc(fss_badengopen)
              call lib$signal (fss_badengopen,%val(1),%val(ios))
          else
              call ct_query_arcv ( , LU_eng, rse, ct_status)

              if ( ct_Status(1) .ne. ctp_normal) then
                   status = %loc(fss_BADQUERY)
                   call lib$signal (fss_BADQUERY, %val(1), %val(CT_STATUS(1)))
              end if
          endif

      elseif (status .and. (rse(1) .eq. ' ') ) then

          open ( unit     = LU_eng,
     .           file     = engfile,
     .           iostat   = ios,
     .           status   = 'OLD',
     .           useropen = CT_CONNECT_READ)
 
          if (ios .ne. 0) then
              status = %loc(fss_badengopen)
              call lib$signal (fss_badengopen,%val(1),%val(ios))
          end if
      endif
c
c Get the engineering records and save the science times
c
      call ct_gmt_to_binary(jstart, start)
      call ct_gmt_to_binary(jstop,  stop)

      do while (status .and. (ct_status(1) .ne. ctp_endoffile) .and.
     .         (num_eng(1) .lt. fss_max_buff)                  .and.
     .         (num_eng(2) .lt. fss_max_buff)                  .and.
     .         (num_eng(3) .lt. fss_max_buff)                  .and.
     .         (num_eng(4) .lt. fss_max_buff) )


         if (rse(1) .ne. ' ') then

             call ct_query_get(, lu_eng, eng_rec, ct_status)
             if ( (ct_status(1) .ne. ctp_normal) .and.
     .            (ct_status(1) .ne. ctp_endoffile)) then
                   status = %loc(fss_badqget)    
                   call lib$signal(fss_badqget, %val(1), %val(ct_status(1)))
             endif

         else

             call ct_read_arcv(, lu_eng, eng_rec, ct_status)
             if ( (ct_status(1) .ne. ctp_normal) .and.
     .            (ct_status(1) .ne. ctp_endoffile)) then
                   status = %loc(fss_badseqrd)    
                   call lib$signal(fss_badseqrd, %val(1), %val(ct_status(1)))
             endif

         endif


         if (status .and. (ct_status(1) .ne. ctp_endoffile) ) then

c
c Get science times from eng records. These will be used for keyed reads.
C Since not all channels may have a time filled in for each eng record, 
c keep number of eng records by channel.
c Number of eng records is only incremented if:
c   1. Time can be converted to an ascii format time via sys$asctim.
c   2. Time is >= start time.
c   3. Time is <= end   time.
c   4. Quality is within limits.
c   5. XCAL position is out.
c.  6. Science mode is in list.
c.  7. MTM mode is in list.
c   8. FAKEIT bit is correct.
c   9. Dihedral temp is low enough.
c

c
c Arrays of temp controller temps are kept because they will 
c be filled into the short science records.
c
             temp_stat=fut_temp_list(eng_rec.en_analog,rawsigs,grtwt,grttrans,
     .                               combswt,singlifg,outtemps,outsigs)
c
c Report if any GRTs had all switches off. Temp_stat = zero if all OK.
c
             if ((temp_stat .gt. 0) .and. (temp_stat .le. 10)) then
                 write(6,'(x,a,i2,a)') 'GRT ',temp_stat,' had all switches off.'
                 write(6,'(x,4a)') 'GRT is ',grt_name(temp_stat),' Time is ',
     &                              eng_rec.ct_head.gmt
             endif

             xcal_pos1 = eng_rec.en_xcal.pos(1)
             xcal_pos2 = eng_rec.en_xcal.pos(2)
             ical_temp = outtemps(2)

             in_bin    = .false.
             do j = 1,num_bins
                if (abs(ical_temp - bins(j)) .le. temptol) in_bin = .true.
             enddo

             i = start_chan
             do while (status .and. (i .le. stop_chan))

                all_eng(i) = all_eng(i) + 1

                adt(1)     = eng_rec.en_head.sci_time(i).bin_time(1)
                adt(2)     = eng_rec.en_head.sci_time(i).bin_time(2)
                sci_mode   = eng_rec.chan(i).up_sci_mode
                speed      = eng_rec.chan(i).xmit_mtm_speed
                length     = eng_rec.chan(i).xmit_mtm_len
                mtm        = (2 * length) + speed 
                iqual      = eng_rec.en_head.dq_summary_flag(i)
                aqual      = eng_rec.en_head.att_summary_flag(i)
                eng_fakeit = eng_rec.chan(i).fakeit
                

                if (check_modes) then 
                    mode_in_list = .false.
                    do j = 1,num_modes
                       if (sci_mode .eq. mpmodes(j)) mode_in_list=.true.
                    enddo
                endif

                if (check_mtm) then 
                    if (mtmlist(mtm+1) .eq. 1) then
                         mtm_in_list = .true.
                    else
                         mtm_in_list = .false.
                    endif
                endif

                if (time_ge(adt, start) .and. time_le(adt,stop)) then
                    in_time_range = .true.
                else
                    in_time_range = .false.
                endif

                if ((iqual .gt. instr_qual).or.(aqual .gt. attit_qual)) then
                    out_qual(i) = out_qual(i) + 1
                elseif ((adt(1) .eq. 0) .and. (adt(2) .eq. 0))          then
                    out_notime(i) = out_notime(i) + 1
                elseif ((.not. in_time_range) .and. (.not. eng_spec) )  then
                    out_time(i) = out_time(i) + 1
                elseif ((xcal_pos1 .ne. fac_xcalout)                    .or.
     .                  (xcal_pos2 .ne. fac_xcalout))                   then
                    out_xcal(i) = out_xcal(i) + 1
                elseif (eng_fakeit .ne. fakeit)                         then
                    out_fake(i) = out_fake(i) + 1
                elseif (outtemps(5) .gt. max_dihed)                     then
                    out_dihed(i) = out_dihed(i) + 1
                elseif (.not. mode_in_list)                             then
                    out_mode(i) = out_mode(i) + 1
                elseif (.not. mtm_in_list)                              then
                    out_mtm(i) = out_mtm(i) + 1
                elseif (.not. in_bin)                                   then
                    out_bin(i) = out_bin(i) + 1
                endif


                if ((in_time_range .or. eng_spec)                         .and.
     .              (xcal_pos1 .eq. fac_xcalout)                          .and.
     .              (xcal_pos2 .eq. fac_xcalout)                          .and.
     .              mode_in_list .and. mtm_in_list .and. in_bin           .and.
     .              time_ge(adt,earliest)                                 .and.
     .              (iqual .le. instr_qual) .and. (aqual .le. attit_qual) .and.
     .              (eng_fakeit .eq. fakeit)                              .and. 
     .              (outtemps(5) .le. max_dihed))                          then

                    num_eng(i) = num_eng(i) + 1
                    if (num_eng(i) .gt. fss_max_buff) then
                         call lib$signal(fss_maxdata)
                         status = %loc(fss_maxdata)
                    endif

                    if (status) then 
                        sci_times(1,i,num_eng(i)) = adt(1)
                        sci_times(2,i,num_eng(i)) = adt(2)

                        xcal(i,num_eng(i)) = outtemps(1)
                        ical(i,num_eng(i)) = outtemps(2)
                        skyh(i,num_eng(i)) = outtemps(3)
                        refh(i,num_eng(i)) = outtemps(4)
                        dihd(i,num_eng(i)) = outtemps(5)
                        bolo(i,num_eng(i)) = outtemps(6+i)
                    endif

                 endif  ! (status and ct_status(1) eq ctp_normal)

                 i = i + chan_int

               enddo    ! (while status and (i le stop_chan))
         endif          ! (status and ct_status(1) .ne. ctp_endoffile)
      enddo             ! (while status and buf_counter le fss_max_buff)
c
c REPORT HOW MANY ENG RECS WERE REJECTED
c
      if (status) then

          j = start_chan
          k = stop_chan
          l = chan_int
          if (report) then 
              stop_lu = lu_rep
              lu_int  = lu_rep - 6
          else
              stop_lu = 6
              lu_int  = 1
          endif

          do lu = 6, stop_lu, lu_int
             write(lu,*) ' '
	     write (lu,fmt='(x,a,t24,7x,4(a2,9x))') 
     .					'CHANNEL ',	(chan_id(i), i=j,k,l)
             write(lu,10) 'INPUT ENG RECS',		(all_eng(i),    i=j,k,l)
             write(lu,10) 'Eng accepted    ',		(num_eng(i),i=j,k,l)
             write(lu,fmt='(x,/,x,a)')  'Eng failure reasons'
             write(lu,10) 'ZERO sci time   ',		(out_notime(i), i=j,k,l)
             write(lu,10) 'Outside time range ',	(out_time(i),   i=j,k,l)
             write(lu,10) 'Quality check ',		(out_qual(i),   i=j,k,l)
             write(lu,10) 'Fakeit ON       ',		(out_fake(i),   i=j,k,l)
             write(lu,10) 'High dihedral   ',		(out_dihed(i),  i=j,k,l)
             write(lu,10) 'MP mode not in list  ',	(out_mode(i),   i=j,k,l)
             write(lu,10) 'MTM mode not in list ',	(out_mtm(i),    i=j,k,l)
             write(lu,10) 'ICAL not in bin ',		(out_bin(i),    i=j,k,l)
             write(lu,10) 'XCAL not out    ',		(out_xcal(i),   i=j,k,l)
             write(lu,*) ' '
10           format(x,a,t24,4(i9,2x))
          enddo
       endif
c
c IF NO ENGINEERING RECORDS FOUND, SIGNAL ERROR
c
      tot_read = 0
      do j = start_chan, stop_chan, chan_int
         tot_read = tot_read + all_eng(j) 
      enddo

      if (tot_read .eq. 0) then
           call lib$signal(fss_noengrecs)
           status = %loc(fss_noengrecs)
      endif
c
c IF ALL ENGINEERING RECORDS REJECTED, SIGNAL ERROR
c
      if (status) then
          tot_good = 0
          do j = start_chan, stop_chan, chan_int
             tot_good = tot_good + num_eng(j) 
          enddo
      
          if (tot_good .eq. 0) then
              call lib$signal(fss_allengrej)
              status = %loc(fss_allengrej)
          endif
      endif
c
c CLOSE ENGINEERING ARCHIVE
c
      if (status .ne. %loc(fss_badengopen)) then
          call CT_CLOSE_ARCV (, LU_eng, ct_status)
          if (ct_status(1) .ne. ctp_normal) then
              status = %loc(fss_badengclose)
              call lib$signal(fss_badengclose, %val(1), %val(ct_status(1)))
          endif
      endif

      fss_read_eng = status

      return
      end
