      integer*4 function fss_close_skymaps(lu_input,lu_output,first_run,write,
     &                                     closeold, jstart, jstop, sky_jstart,
     &                                     in_def, eng_def)

C------------------------------------------------------------------------
C    PURPOSE:  Close the input and output skymaps.
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION:  status = fss_close_skymaps(lu_input,lu_output,first_run,write)
C
C    INPUT PARAMETERS:       I*4   lu_output     The output skymap LU
C                            I*4   lu_input      The input  skymap LU
C                            I*4   first_run     First run flag
C                            I*4   write         Flag to write skymap
C                            I*4   closeold      CLOSEOLD qualifier
C                            C*14  jstart        Jstart from command line
C                            C*14  jstop         Jstop from command line
C                            I*4   sky_jstart(2) Jstart of input skymap
C                            L*4   in_def        Was Input skymap defaulted ?
C                            L*4   eng_def       Was Eng file defaulted ?
C
C    OUTPUT PARAMETERS:      None.  
C 
C    SUBROUTINES CALLED:     csa_close_skymap
C
C    COMMON VARIABLES USED:  None.
C 
C    INCLUDE FILES:          None.
C
C----------------------------------------------------------------------
C                    PDL for FSS_CLOSE_SKYMAPS
c                    D. Bouler, JAN, 1991
c
c  set jstart and jstop times of skymaps
c
c  if (.not. first_run) call csa_close_skymap to close input skymap
c 
c  if (write) call csa_close_skymap to close output skymap 
c
c  fss_close_skymaps = status
c
c  return
c  end (pdl)
c
c----------------------------------------------------------------------
C Changes:
C
C
C----------------------------------------------------------------------

      implicit none
c
c Function declarations
c
      integer*4      csa_close_skymap
      integer*4      csa_set_time_range
      integer*4      lib$free_lun
C
C Variable declarations
C
      integer*2      index_level /6/       !Number of index levels for skymaps

      integer*4      i,j,k,l               !Loop counters
      integer*4      status                !General status return
      integer*4      lu_output             !Output skymap LU
      integer*4      lu_input              !Input  skymap LU
      integer*4      first_run             !Flag for first run of month
      integer*4      write                 !Flag to write skymap
      integer*4      closeold              !CLOSEOLD qualifier
      integer*4      sky_jstart(2)         !Jstart of input skymap
      integer*4      start(2)              !Jstart time of output skymap
      integer*4      stop(2)               !Jstop  time of output skymap

      character*14   jstart                !Jstart from command line
      character*14   jstop                 !Jstop  from command line

      logical*4      in_def                !Was input skymap defaulted ?
      logical*4      eng_def               !Was Eng file defaulted ?
c
c Declare externals
c
      external       fss_normal
c*
c******  Begin code *******************************************
c*
      fss_close_skymaps = %loc(fss_normal)         
      status            = %loc(fss_normal)         
c*
c* SET JSTART AND JSTOP TIMES OF OUTPUT SKYMAP
c*
      if ((first_run .or. closeold) .and. write) then
          call ct_gmt_to_binary(jstart, start)         
          call ct_gmt_to_binary(jstop , stop)         
          status = csa_set_time_range(lu_output, start, stop)
          if (.not. status) call lib$signal(%val(status))
      elseif ((.not. first_run) .and. (.not. closeold)  .and.
     &         in_def .and. eng_def .and. write)         then
          call ct_gmt_to_binary(jstop , stop)         
          status = csa_set_time_range(lu_output, sky_jstart, stop)
          if (.not. status) call lib$signal(%val(status))
      endif
c*
c* CLOSE INPUT SKYMAP IF NOT FIRST RUN 
c*
       if ( .not. first_run ) then
           status = csa_close_skymap( lu_input, index_level)
           if (.not. status) call lib$signal(status)
           status = lib$free_lun(lu_input)
           if (.not. status) call lib$signal(status)
       end if
c*
c* CLOSE OUTPUT SKYMAP IF WRITE HAS BEEN DONE
c*
      if (write) then
          status = csa_close_skymap( lu_output, index_level)
          if (.not. status) call lib$signal(status)
          status = lib$free_lun(lu_output)
          if (.not. status) call lib$signal(status)
      endif

      fss_close_skymaps = status

      return
      end
