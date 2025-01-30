
        integer*4 function fss_open_skymaps(lu_input, lu_output, first_run, 
     .                                      closeold, oldest, 
     .                                      in_ext,   out_ext,   write,
     .                                      jstart, jstop, chan, fdq_ext,
     .                                      sky_jstart,in_def,out_def,eng_ext)

C------------------------------------------------------------------------
C    PURPOSE: Opens the input and output skymaps.
C
C    AUTHOR: D. Bouler, STX, Jan , 1990
C
C    INVOCATION:  status = fss_open_skymaps(lu_input, lu_output, first_run, 
C                                           closeold, oldest, 
C                                           in_ext,   out_ext,   write,
C                                           jstart, jstop, chan, fdq_ext,
C                                           sky_jstart,in_def,out_def,eng_ext)
C
C    INPUT PARAMETERS:    i*4   write     whether to write output skymap
C                         i*4   first_run whether input skymap is used
C                         i*4   closeold  Include uncoadded groups from skymap?
C                         c*14  oldest    earliest acceptable time from skymap
C                         c*128 in_ext    input skymap extension
C                         c*128 out_ext   output skymap ext 
C                         c*14  jstart    start process time
C                         c*14  jstart    start process time
C                         I*4   chan      channel being processed
C                         c*3   fdq_ext   FDQ file extension ('ED_' for edit)
C                         L*4   in_def    Was input skymap defaulted?
C                         L*4   out_def   Was output skymap defaulted?
C
C    OUTPUT PARAMETERS:   I*4   lu_input      LU for input skymap
C                         I*4   lu_output     LU for output skymap
C                         I*4   sky_jstart(2) JSTART of input skymap
C                         c*128 eng_ext       Eng file extension
C
C    SUBROUTINES CALLED:    cct_query_ttg_catalog
C                           str$trim
C                           fut_clean_catalog
c			    cct_query_catalog
c			    ct_binary_to_gmt
c			    lib$signal
c
C    COMMON VARIABLES USED: None.
C 
C    INCLUDE FILES:         fut_params
C                           $ssdef
C                           ct$library:ctuser.inc
C                           fss_include
C                           CCT_Query_TTG_Catalog_Record
c			    cct_query_catalog_record.txt
C
C
C----------------------------------------------------------------------
C                      PDL for FSS_open_skymaps
c                      D. BOULER, Jan , 1990
c
c  Build eng file extension 
c
c  if (.not. first_run) then
c       if (input skymap extension specified) then
c           build input skymap name from extension
c       elseif (input skymap name defaulted) then
c          call cct_query_ttg_catalog to get latest input skymap name
c      endif
c
c      status = lib$get_lun(lu_input)
c      if (.not. status) then
c           signal error
c      else
c           call csa_open_skymap to open input skymap
c           if (.not. status) signal error
c      endif          
c  endif          
c
c  if (write) then
c      if (output skymap defaulted) then
c	 if ( first_run ) then
c          build default skymap name using current JSTART and JSTOP
c	 else
c          build default skymap name using input skymap JSTART and current JSTOP
c	 endif
c      else
c          build output skymap name using extension
c      endif
c
c      status = lib$get_lun(lu_output)
c      if (.not. status) then
c           signal error
c      else
c           call csa_open_skymap to open output skymap
c           if (.not. status) signal error
c      endif
c  endif
c
c  fss_open_skymaps = status
c
c  return
c  end (pdl)
c
c-------------------------------------------------------------------------------
C Changes:
C 
C     SPR 9962, Tilak Hewagama, Hughes STX 10-Sept-1992
c          Use the start time from the input skymap catalog and the stop time 
c          from the current processing session to build the default output 
c          skymap file name.
C
C
C-------------------------------------------------------------------------------

      implicit none

      Include '($SSDef)'
      Include 'CT$LIBRARY:CTUser.Inc'
      Include '(CCT_Query_TTG_Catalog_Record)'
      Include '(CCT_Query_Catalog_Record)'
      Include '(FUT_Params)'
      Include '(fss_include)'
                                             
      Dictionary 'CCM_CME_Catalog_Entry'

      Record /CCM_CME_CATALOG_ENTRY/ cats(50)
      Record /QUERY_TTG_CATALOG/     query
      Record /QUERY_CATALOG/         query_cat
c*
c* FUNCTION DECLARATIONS
c*
      Integer*4      CCT_Query_TTG_Catalog
      Integer*4      CCT_Query_Catalog
      integer*4      lib$get_lun
      integer*4      cli$present
      integer*4      cli$get_value
      logical*4      time_lt
      logical*4      time_gt
      logical*4      time_le
      logical*4      time_ge
      integer*4      fut_clean_catalog
C*
C* VARIABLE DECLARATIONS
C*
      Integer*2      cats_num
      Integer*2	     ct_status(20)

      integer*4      i,j,k,l
      integer*4      status
      integer*4      cstatus
      integer*4      lu_sci
      integer*4      lu_input
      integer*4      lu_output
      integer*4      first_run
      integer*4      closeold
      integer*4      chan
      integer*4      ios
      integer*4      start(2)              ! Jstart in binary format
      integer*4      stop(2)               ! jstop  in binary format
      integer*4      earliest(2)           ! oldest in binary format
      integer*4      inp_astart(2)         ! input skymap start in binary format
      integer*4      activity(2)           ! Activity time 
      integer*4      inp_len               ! Length of input  file name
      integer*4      out_len               ! Length of output file name
      integer*4      fblen1                ! filename buffer length
      integer*4      fblen2                ! filename buffer length
      integer*4      iday                  ! Day extracted from jstart
      integer*4      iyear                 ! Year extracted from jstart
      integer*4      write                 ! Whether to write output skymap
      integer*4      last                  ! Number of skymap from last run 
      integer*4      period                ! Position of "." in skymap name
      integer*4      purge     /0/         ! Flag used in FUT_CLEAN_CATALOG;
                                           ! purge duplicates if = 1
      integer*4      last_map              ! Number of most recent skymap
      integer*4      sky_jstart(2)         ! Jstart of input skymap
      integer*4      eng_len               ! Length of eng extension

      character*2    achan                 ! Ascii channel
      character*2    cyear                 ! Ascii year of start time

      character*3    cday                  ! Day extracted from jstart
      character*3    fdq_ext               ! FDQ file extension (ED_ for edit)

      character*14   mstart                ! Start of month
      character*14   jstart
      character*14   jstop
      character*14   inp_jstart
      character*14   oldest                ! Earliest acceptable time
      character*14   t_accept              ! Earliest acceptable time

      character*64    filebuff1            ! Filename buffer
      character*64    filebuff2            ! Filename buffer

      character*128 in_ext                 ! Input skymap  extension
      character*128 out_ext                ! Output skymap extension
      character*128 eng_ext                ! Eng file      extension
      character*128 inskymap               ! Full input skymap name 
      character*128 outskymap              ! Full output skymap name

      logical*4     in_def                 ! Was input  skymap defaulted ?
      logical*4     out_def                ! Was output skymap defaulted ?

      external       fss_ctquerycat
      external       fss_nocatrecs
      external       fss_normal
      external       fss_noinfile
      external       csa_open_skymap
      external       fss_dupskymap
      external       fss_allmapsrej1
      external       fss_allmapsrej2
      external       cct_normal
      external       cct_q_no_cat_entry
c*
c******  Begin code ********
c*
       fss_open_skymaps = %loc(fss_normal)         
       status           = %loc(fss_normal)
       achan            = fac_channel_ids(chan)
       sky_jstart(1)    = 0
       sky_jstart(2)    = 0
c*
c* Build eng file extension.
c*
       call str$trim(eng_ext, eng_ext, eng_len)

       if (eng_ext .eq. ' ') then
           eng_ext = fdq_ext // jstart(:11) // ';' // jstop(:11)
       else
           period = index(eng_ext, '.')
           if (period .gt. 0) then
               eng_ext = eng_ext(period+1:eng_len)
           else
               eng_ext = eng_ext(:eng_len)
           endif
       endif
c*
c* Get input skymap name if defaulted. Go back 30 days to search catalog.
c*
       if (status .and. (.not. first_run)) then

           if  (in_def) then

                in_ext = ' '

                cyear = jstart(1:2)
                cday  = jstart(3:5)

                read(cyear,fmt='(i2)',iostat=ios) iyear
                if (ios .ne. 0) then
                    call errsns(ios,,,,status)
                    call lib$signal(%val(status))
                endif             

                read(cday,fmt='(i3)',iostat=ios) iday
                if (ios .ne. 0) then
                    call errsns(ios,,,,status)
                    call lib$signal(%val(status))
                endif             

                if (status .and. (iday .lt. 31)) then              
                    iyear = iyear - 1
                    iday  = 365 + iday - 30 
                elseif (status .and. (iday .ge. 31)) then
                    iday = iday - 30
                endif

                if (status) then
                    write(cday,fmt='(i3.3)',iostat=ios) iday
                    if (ios .ne. 0) then 
                        call errsns(ios,,,,status)
                        call lib$signal(%val(status))
                    endif

                    write(cyear,fmt='(i2.2)',iostat=ios) iyear
                    if (ios .ne. 0) then 
                        call errsns(ios,,,,status)
                        call lib$signal(%val(status))
                    endif
                 endif

                 if (status) then
                     mstart = cyear(1:2) // cday(1:3) // jstart(6:14)

                     call ct_gmt_to_binary(mstart, query.start_time)
                     call ct_gmt_to_binary(jstart, query.stop_time)

                     query.archive_id    = 'CSDR$FIRAS_IN'
                     query.dataset_name  = 'FSS_SSSKY_' // achan

                     status = CCT_Query_TTG_Catalog(query, cats, cats_num)
                 endif

                 If (.Not. status) then
                     status = %loc(fss_ctquerycat)
                     call lib$signal(fss_ctquerycat, %val(1), %val(status))
                 elseif (cats_num .eq. 0) then 
                     status = %loc(fss_nocatrecs)
                     call lib$signal(fss_nocatrecs, %val(1),query.dataset_name)
                 else
                     status = fut_clean_catalog(cats, cats_num, purge)

                     if (.not. status) then
                          call lib$signal(%val(status))
                     else
                         last_map    = -1
                         activity(1) = 0
                         activity(2) = 0
                         do i = 1, cats_num
                            if (time_gt(cats(i).activity_time, activity)   .and.
     &                          time_le(cats(i).final_time,query.stop_time))then
                                activity(1) = cats(i).activity_time(1)
                                activity(2) = cats(i).activity_time(2)
                                last_map    = i
                            endif
                         enddo

                         if (last_map .eq. -1) then
                             call lib$signal(fss_allmapsrej1)
                             call lib$signal(fss_allmapsrej2)
                             status = %loc(fss_allmapsrej1)
                         else
                             sky_jstart(1) = cats(last_map).initial_time(1)
                             sky_jstart(2) = cats(last_map).initial_time(2)
                        
                             if (cats(last_map).filename_extension .eq.
     .                           cats(last_map-1).filename_extension)   then
                                 call lib$signal(fss_dupskymap)
                             endif
                         endif
                     endif

                     in_ext = cats(last_map).filename_extension
                     call str$trim(in_ext, in_ext, inp_len)
                     inskymap = 'CSDR$FIRAS_IN:FSS_SSSKY_' // achan // '.' // 
     .                           in_ext(:inp_len)
                 endif

           elseif (.not. in_def) then

                 call str$trim(in_ext, in_ext, inp_len)
                 period = index(in_ext, '.')
                 if (period .gt. 0) then
                     inskymap = 'CSDR$FIRAS_IN:FSS_SSSKY_' // achan // '.' // 
     .                           in_ext(period+1:inp_len)
                     in_ext   =  in_ext(period+1:inp_len)
                 else
                     inskymap = 'CSDR$FIRAS_IN:FSS_SSSKY_' // achan // '.' // 
     .                           in_ext(:inp_len)
                     in_ext   =  in_ext(:inp_len)
                 endif

           endif   !  (in_def)

       endif       !  (status .and. (.not. first_run))
c*
c* OPEN THE INPUT SKYMAP.
c*
       if ( status .and. (.not. First_Run) ) then
       
            status = lib$get_lun(LU_input)
            if (.not. status) call lib$signal(%val(status))
   
            call str$trim(inskymap, inskymap, inp_len)

            if (status) then
                open ( unit       = LU_input,
     .                 file       = inskymap(1:inp_len),
     .                 iostat     = ios,
     .                 status     = 'OLD',
     .                 form       = 'unformatted',
     .                 recordtype = 'fixed',
     .                 readonly,
     .                 useropen   = csa_open_skymap)

                 if (ios .ne. 0) then
                     call errsns(ios,,,,status)
                     call lib$signal(%val(status))
                 endif
            end if
        end if
c*
c* Build output skymap name if defaulted; use JSTART/JSTOP
c* if an input skymap is not used, else use JSTART from input
c* skymap catalog and JSTOP from current session. If the file
c* name is not defaulted, then look for "." in file name. 
c* Accept everything after the "." if there is one; otherwise, 
c* just use as an extension.
c*
       if (write) then

           if (out_def) then
	     if ( first_run ) then
c            Use current session JSTART and JSTOP for output file name
               outskymap = 'CSDR$FIRAS_OUT:FSS_SSSKY_' // achan // '.' // 
     .                      fdq_ext // jstart(:7) // '_' // jstop(:7)
               out_ext   =  fdq_ext // jstart(:7) // '_' // jstop(:7)
	     else if ( ( .not. first_run ) .and. closeold ) then
c            Use OLDEST as JSTART
               call ct_gmt_to_binary(oldest,earliest)
               call ct_binary_to_gmt(earliest,t_accept)
               outskymap = 'CSDR$FIRAS_OUT:FSS_SSSKY_' // achan // '.' // 
     .                      fdq_ext // t_accept(:7) // '_' // jstop(:7)
               out_ext   =  fdq_ext // t_accept(:7) // '_' // jstop(:7)
	     else if ( ( .not. first_run ) .and. ( .not. closeold ) ) then
c            Define the archive name and query the Cobetrieve catalog.

		query_cat.archive_id = 'CSDR$FIRAS_IN'
		query_cat.filename   = inskymap(15:inp_len)

		cstatus = cct_query_catalog (query_cat, cats(1))
		call str$trim (filebuff1, query_cat.archive_id, fblen1)
		call str$trim (filebuff2, query_cat.filename, fblen2)

		if (cstatus .eq. %loc(cct_q_no_cat_entry)) then
		   fss_open_skymaps = %loc(fss_nocatrecs)
		   call lib$signal (fss_nocatrecs, %val(2), filebuff1(1:fblen1),
     .							 filebuff2(1:fblen2))
		   return
		elseif (cstatus .ne. %loc(cct_normal)) then
		   fss_open_skymaps = %loc(fss_ctquerycat)
		   call lib$signal (fss_ctquerycat,%val(3), filebuff1(1:fblen1),
     .					 filebuff2(1:fblen2), %val(cstatus))
		   return
		endif

c            If OK, then get the start time for the skymap catalog and form name

		inp_astart(1) = cats(1).initial_time(1)
		inp_astart(2) = cats(1).initial_time(2)
		call ct_binary_to_gmt (inp_astart,inp_jstart)

		outskymap = 'CSDR$FIRAS_OUT:FSS_SSSKY_' // achan // '.' // 
     .                      fdq_ext // inp_jstart(:7) // '_' // jstop(:7)
		out_ext   =  fdq_ext // inp_jstart(:7) // '_' // jstop(:7)

	     endif                     ! first_run default file name
           else                   ! End branch for default name
               call str$trim(out_ext, out_ext, out_len)
               period = index(out_ext, '.')
               if (period .gt. 0) then
                   outskymap = 'CSDR$FIRAS_OUT:FSS_SSSKY_' // achan // '.' //
     .                          out_ext(period+1:out_len)
                   out_ext    = out_ext(period+1:out_len)
               else
                   outskymap = 'CSDR$FIRAS_OUT:FSS_SSSKY_' // achan // '.' //
     .                          out_ext(:out_len)
                   out_ext   =  out_ext(:out_len)
               endif
           endif                 ! (out_def) End branch for user specified name
           call str$trim(outskymap, outskymap, out_len)
        endif
c*
c* OPEN OUTPUT SKYMAP FILE
c*
        if (status .and. write) then

            status = lib$get_lun(LU_Output)
            if (.not. status) call lib$signal(%val(status))

            if (status) then
                open ( unit       = LU_Output,
     .                 file       = outskymap(1:out_len),
     .                 iostat     = ios,
     .                 status     = 'NEW',
     .                 recl       = fss_lw_recsize,
     .                 form       = 'unformatted',
     .                 recordtype = 'fixed',
     .                 useropen   = csa_open_skymap)

                if (ios .ne. 0) then 
                    call errsns(ios,,,,status)
                    call lib$signal(%val(status))
                endif
            endif
        endif

        fss_open_skymaps = status

        return
        end
