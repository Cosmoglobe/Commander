      integer*4 function fss_read_config(minimum, grtwt, grttrans)

C------------------------------------------------------------------------
C    PURPOSE: This routine reads the configuration files FEX_MINCOADD,
C             FEX_GRTRAWWT, and FEX_GRTTRANS. The minimum ifgs in a group is 
C             from FEX_MINCOADD, the flags for which weights are used for GRTs
C             comes from FEX_GRTRAWWT, and the GRT transition temps come from 
C             FEX_GRTTRANS.
C             
C    AUTHOR: D. Bouler, STX, Dec, 1990
C
C    INVOCATION:  status = fss_read_config(minimum, grtwt, grttrans)
C
C    INPUT PARAMETERS:         None.
C
C    OUTPUT PARAMETERS:        I*4    minimum(4)   Min ifgs in coadd
C                              Record grtwt        Grt weights
C                              Record grttrans     Grt transition temps
C
C    SUBROUTINES CALLED:       cct_open_config
C                              cct_get_config_tod
C                              cct_close_config
C
C    COMMON VARIABLES USED:    None.
C 
C    INCLUDE FILES:            ct$library:ctuser.inc
C                              cct_get_config
C
C----------------------------------------------------------------------
c
c                  PDL for FSS_READ_CONFIG
c                  D. Bouler, DEC, 1990
c
c   Set return status to fss_normal
c
c   Open configuration files 
c   if (.not. status) signal error
c
c   if (status) then
c       Read minimum number of coadds in a group
c       if (.not. status) signal error
c   endif
c   
c   if (status) then
c       Read grt weight information
c       if (.not. status) signal error
c   endif
c   
c   if (status) then
c       Read grt transition temp information
c       if (.not. status) signal error
c   endif
c   
c   close configuration files
c   if (.not. status) signal error
c
c   fss_read_config = status
c
c   return
c   end (pdl)
c*
c*-----------------------------------------------------------------------------
c*    Changes:
c*
c*
c*-----------------------------------------------------------------------------

      implicit none
c*
c*    Include Files:
c*
      include 'CT$LIBRARY:CTUSER.INC'      !COBETRIEVE definitions
      include '(cct_get_config)'
c*
c*    Functions
c*
      integer*4 CCT_Get_Config_TOD
      integer*4 CCT_Open_Config
      integer*4 CCT_Close_Config
c*
c*    Externals
c*
      external fss_normal                  !Everything OK signal
c*
c*    Variables:
c*
      integer*2 CT_STATUS(20)              !CT RETURN STATUS

      integer*4 i,j,k,l                    !Loop counters
      integer*4 minimum(4)                 !Min number of ifgs in a coadd
      integer*4 ios                        !Input/Output status
      integer*4 status                     !general status
      integer*4 Lu_Config(3)               !logical unit number
      integer*4 index(3)                   !initial cache pointers
      integer*4 ncache                     !number of caches to use
      integer*4 Ref_Count(3)               !number of config datasets
      integer*4 ndsets                     !number of config datasets
      integer*4 size(3)                    !size of config datasets
      integer*4 start(2)                   !Start time for config (870010000)
      integer*4 stop(2)                    !Stop  time for config (993650000)

      logical*1 New_Segment(3)             !flag for new config segments
      
      character*1 Access_Mode/' '/         !dataset access mode

      character*14 jstart                  !GMT start time
      character*14 jstop                   !GMT stop  time

      character*32 name(3)                 !name of the datasets
c*
c* Declare dictionary and record structures
C*
      dictionary   'fex_mincoadd'
      dictionary   'fex_grtrawwt'
      dictionary   'fex_grttrans'

      record /config_status/  stat(3)

      structure  /config_rec/
        record /fex_mincoadd/  fex_mincoadd
        record /fex_grtrawwt/  fex_grtrawwt
        record /fex_grttrans/  fex_grttrans
      end structure

      record /config_rec/   config_rec
      record /fex_grtrawwt/ grtwt
      record /fex_grttrans/ grttrans

c*
c***  BEGIN CODE  ******************************************************
c*
      fss_read_config = %loc(fss_Normal)
      status  = %loc(fss_Normal)
      name(1) = 'csdr$firas_ref:fex_mincoadd'
      name(2) = 'csdr$firas_ref:fex_grtrawwt'
      name(3) = 'csdr$firas_ref:fex_grttrans'
      ncache  = 3
      ndsets  = 3
      size(1) = 128
      size(2) = 256
      size(3) = 384
c*
c* OPEN CONFIG FILE 
c*
      jstart = '88320000000000'
      jstop  = '99360000000000'
      call ct_gmt_to_binary(jstart,start)
      call ct_gmt_to_binary(jstop,stop)

      if (status) then
          status = cct_open_config(Start, Stop, ndsets, name, size,
     .                             Access_Mode, Ncache, Lu_Config,
     .                             Index, Stat, Ref_count)
          if (.not. status) call lib$signal(%val(status))
      endif
c*
c* READ RUNTIME PARAMETERS FROM CONFIG FILE 
c*
      if (status) then
          status = cct_get_config_tod(start, ndsets, size,
     .                                lu_config, index, config_rec, 
     .                                new_segment, stat)
          if (.not. status) then
              call lib$signal(%val(status))
          else
              do i = 1,4
                 minimum(i) = config_rec.fex_mincoadd.min_ifg_coadd(i)
              enddo
              grtwt      = config_rec.fex_grtrawwt
              grttrans   = config_rec.fex_grttrans
          endif
      endif
c*
c* Close configuration files
c*
      if (status) then
          status = CCT_Close_Config(Ndsets, LU_Config, Index)
          if (.not. status) call lib$signal(%val(status))
      endif

      fss_read_config = status

      return
      end
