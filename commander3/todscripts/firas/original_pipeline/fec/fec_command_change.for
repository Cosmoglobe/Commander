      integer*4 function FEC_Command_Change (
     R                                       Index,
     R                                       Array)
c
c     By Reid Wilson, SASC Technologies Inc.,  13-AUG-1986
c
c     Description:
c         FEC_Command_Change checks whether any of the side A and B
c         temperatures have been commanded to change.  A true value
c         is returned if a change occurs, otherwise false is returned.
c
c     Format of Call:
c         Ret_Code = FEC_Command_Change (
c    R                                   Index,
c    R                                   Array)
c
c     Input Parameters:
c         Index            = index of engineering record in Array to 
c                            check for command change.
c                            integer*4
c         Array            = array of engineering records
c                            RECORD /FDQ_ENG_DATA/ Array(*)
c
c     Output Parameters:
c         Ret_Code         = TRUE if a change command occurred, otherwise FALSE
c                            integer*4
c
c     Modules Called:
c         None
c
c     Include Files: 
c         FEC.PAR          = parameters used by FEC (FEC_TRUE, FEC_FALSE)
c         FIRAS_REC_DEF    = definition of engineering record
c
c     ******* BEGIN *******
c
      implicit none
c
c     Include files:
c
      include '(FEC_INC)'             ! contains TRUE/FALSE values to return
c
c     Formal Parameters:
c
      integer*4 Index                 ! index to element of Array to check    
      integer*4 I
      integer*4 Flag
c
      DICTIONARY 'FDQ_ENG'
      RECORD /FDQ_ENG/ Array(*)       ! array of engineering records
c
c     *** BEGIN ***
c
      Flag = 0
      I = Index
c     if index = 1 then cannot check previous record: return FALSE
      If (I .eq. 1) then
       FEC_Command_Change = FEC_FALSE
       return
      Endif

c     If index gt 1, then look for command changes between current and previous
c        records. If any change, flag the change and return.

      If (I .gt. 1) then

         if(Array(I).EN_STAT.BOL_CMD_BIAS(1) .ne.
     +      Array(I-1).EN_STAT.BOL_CMD_BIAS(1) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if(Array(I).EN_STAT.BOL_CMD_BIAS(2) .ne.
     +      Array(I-1).EN_STAT.BOL_CMD_BIAS(2) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if(Array(I).EN_STAT.BOL_CMD_BIAS(3) .ne.
     +      Array(I-1).EN_STAT.BOL_CMD_BIAS(3) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if(Array(I).EN_STAT.BOL_CMD_BIAS(4) .ne.
     +      Array(I-1).EN_STAT.BOL_CMD_BIAS(4) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_1,12,2) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_1,12,2) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_4,0,16) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_4,0,16) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_5,12,2) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_5,12,2) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_8,0,16) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_8,0,16) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_9,12,2) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_9,12,2) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_12,0,16) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_12,0,16) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_13,12,2) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_13,12,2) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.STAT_WORD_16,0,16) .ne.
     +      IBITS(Array(I-1).EN_STAT.STAT_WORD_16,0,16) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.INT_REF_TEMP_A,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.INT_REF_TEMP_A,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.INT_REF_TEMP_B,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.INT_REF_TEMP_B,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.REF_HRN_TEMP_A,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.REF_HRN_TEMP_A,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE 
               return
         endif
         if( IBITS(Array(I).EN_STAT.REF_HRN_TEMP_B,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.REF_HRN_TEMP_B,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.SKY_HRN_TEMP_A,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.SKY_HRN_TEMP_A,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.SKY_HRN_TEMP_B,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.SKY_HRN_TEMP_B,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.EXT_CAL_TEMP_A,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.EXT_CAL_TEMP_A,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif
         if( IBITS(Array(I).EN_STAT.EXT_CAL_TEMP_B,0,12) .ne.
     +      IBITS(Array(I-1).EN_STAT.EXT_CAL_TEMP_B,0,12) ) then
               Flag = 1
               FEC_Command_Change = FEC_TRUE
               return
         endif

      Endif
c
c     if the current values are different than those in the previous record, 
c     then return TRUE
      IF ( Flag .eq. 1) Then
          FEC_Command_Change = FEC_TRUE
      else
c     no change: return FALSE
          FEC_Command_Change = FEC_FALSE
      End if
c
      return
      end !FEC_Command_Change
