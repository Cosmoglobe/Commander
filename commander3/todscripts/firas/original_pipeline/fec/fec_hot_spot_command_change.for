      integer*4 function FEC_Hot_Spot_Command_Change (side,hot_start,hot_end,
     R                                       Index,
     R                                       Array)
c
c     By H. Wang, STX,  22-May-1990
c
c     Description:
c         FEC_Hot_Spot_Command_Change checks whether any of the side A and B
c         Hot_Spot temperature have been commanded to change.  A true value
c         is returned if a change occurs, otherwise false is returned.
c
c     Format of Call:
c         Ret_Code = FEC_Hot_Spot_Command_Change (side,hot_start,hot_end,
c    R                                   Index,
c    R                                   Array)
c
c     Input Parameters:
c         Side             = 1  hot spot A side
c                            2  hot spot B side
c         Index            = index of engineering record in Array to 
c                            check for command change.
c                            integer*4
c         Array            = array of engineering records
c                            RECORD /FDQ_ENG_DATA/ Array(*)
c
c     Output Parameters:
c         hot_start(2)     = the flag to indicate hot spot event is start     
c         hot_end(2)       = the flag to indicate hot spot event is end     
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
c     Changes:
c
c     PDL:
c      if hot spot not starting    
c      then 
c        if the fdq_eng record the hot spot temperature is not
c          zero
c        then
c           set fec_hot_spot_command_change to true
c           set hot_spot_starting to true.
c        endif
c      else if( hot spot starting)
c      then
c        if the fdq_eng record, the hot spot temperature is zero
c        then
c           set fec_hot_spot_command_change to true
c           set hot_spot_ending to .true.
c        endif
c      else
c        set fec_hot_spot_command_change to false
c      endif
c      return
c      end
c 
c
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
      integer*4 side
      integer*4 Index                 ! index to element of Array to check    
      logical*1 hot_start(2)
      logical*1 hot_end(2)
      logical*1 side_start(2)/2*.false./
      logical*1 side_end(2)/2*.false./             
      DICTIONARY 'FDQ_ENG'
      RECORD /FDQ_ENG/ Array(*)       ! array of engineering records
c
c     *** BEGIN ***
c
c     if index = 1 then cannot check previous record: return FALSE
      if (Index .eq. 1) then 
         FEC_Hot_Spot_Command_Change = FEC_FALSE
         return
      endif
      if (Index .gt. 1) then

c         If hot spot commanded OFN: Set hot_start(side) and return TRUE

          if(Array(Index).EN_STAT.HOT_SPOT_CMD(side) .eq. 1
     +       .and. Array(Index-1).EN_STAT.HOT_SPOT_CMD(side) .eq. 0) then
                   hot_start(side) = .true.
                   hot_end(side) = .false.
                   FEC_Hot_Spot_Command_Change = FEC_TRUE
                   return
          endif

c         If hot spot commanded OFF: Set hot_end(side) and return TRUE

          if(Array(Index).EN_STAT.HOT_SPOT_CMD(side) .eq. 0
     +       .and. Array(Index-1).EN_STAT.HOT_SPOT_CMD(side) .eq. 1) then
                   hot_end(side) = .true.
                   hot_start(side) = .false.
                   FEC_Hot_Spot_Command_Change = FEC_TRUE
                   return
          endif

c         If no hot spot command change: return FALSE

          if(Array(Index).EN_STAT.HOT_SPOT_CMD(side) .eq. 
     +       Array(Index-1).EN_STAT.HOT_SPOT_CMD(side)) then
                   FEC_Hot_Spot_Command_Change = FEC_FALSE
                   return
          endif

      End if
c
      return
      end !FEC_Hot_Spot_Command_Change
