function FDQ_IDX_st, Nstruct
 
  common FDQ_IDX_st, defined
 
  if N_elements( defined ) NE 1 then begin
 
       header = { HEADER                           , $
                  GMT_START_TIME: bytarr( 14 )     , $
                  GMT_END_TIME: bytarr( 14 )       , $
                  BNRY_START_TIME: Lonarr( 2 )     , $
                  BNRY_END_TIME: Lonarr( 2 )       , $
                  SC_START_TIME: bytarr( 6 )       , $
                  SC_END_TIME: bytarr( 6 )         , $
                  DATASET_ID: 0                    , $
                  IDXHEAD_SPARES: bytarr( 6 )          }
 
       xcal = { XCAL             , $
                POS: intarr( 2 )     }
 
       ichan = { ICHAN                   , $
                  UP_SCI_MODE: 0B        , $
                  FAKEIT: 0B             , $
                  UP_ADDS_GROUP: 0B      , $
                  UP_SWEEPS_IFG: 0B      , $
                  MTM_SPEED_XMIT: 0B     , $
                  MTM_LENGTH_XMIT: 0B    , $
                  SCIENCE_GAIN: 0        , $
                  BOL_CMD_BIAS: 0B       , $
                  BOL_RDOUT_VOLT: 0B     , $
                  CHAN_SPARE1: 0B        , $
                  CHAN_SPARE2: 0B            }

       grts = { GRTS                   , $
                XCAL: 0                , $
                SKYH: 0                , $
                REFH: 0                , $
                ICAL: 0                , $
                DIHEDRAL_RIGHT: 0      , $
                BOL_ASSEMBLY_RH: 0     , $
                BOL_ASSEMBLY_RL: 0     , $
                BOL_ASSEMBLY_LH: 0     , $
                BOL_ASSEMBLY_LL: 0     , $
                MIRROR_M1_RIGHT: 0     , $
                CAL_RESISTOR_1: 0      , $
                CAL_RESISTOR_2: 0      , $
                CAL_RESISTOR_3: 0      , $
                CAL_RESISTOR_4: 0      , $
                XCAL_SEGMENT: 0        , $
                COLLIM_C1_RIGHT: 0         }
 
       temp_ctrl = { TEMP_CTRL                , $
                     TMP_CTRL_ICAL_A: 0       , $
                     TMP_CTRL_REFH_A: 0       , $
                     TMP_CTRL_SKYH_A: 0       , $
                     TMP_CTRL_XCAL_A: 0       , $
                     TMP_CTRL_ICAL_B: 0       , $
                     TMP_CTRL_REFH_B: 0       , $
                     TMP_CTRL_SKYH_B: 0       , $
                     TMP_CTRL_XCAL_B: 0       , $
                     TMP_CTRL_SPARES: bytarr( 4 ) }
 
       stat_mon = { STAT_MON                     , $
                    MTM_OPER_MODE: 0B            , $
                    ICAL_TEMP_CNTRL: 0B          , $
                    REFH_TEMP_CNTRL: 0B          , $
                    SKYH_TEMP_CNTRL: 0B          , $
                    XCAL_TEMP_CNTRL: 0B          , $
                    DITHER_CHAN_H: 0B            , $
                    DITHER_CHAN_L: 0B            , $
                    ROM_POWER_H: 0B              , $
                    ROM_POWER_L: 0B              , $
                    MTM_LATCH: 0B                , $
                    XCAL_LATCH: 0B               , $
                    ICAL_CURR_RANGE: 0B          , $
                    REFH_CURR_RANGE: 0B          , $
                    SKYH_CURR_RANGE: 0B          , $
                    XCAL_CURR_RANGE: 0B          , $
                    ICAL_INTEG_GAIN: 0B          , $
                    REFH_INTEG_GAIN: 0B          , $
                    SKYH_INTEG_GAIN: 0B          , $
                    XCAL_INTEG_GAIN: 0B          , $
                    ICAL_PROP_GAIN: 0B           , $
                    REFH_PROP_GAIN: 0B           , $
                    SKYH_PROP_GAIN: 0B           , $
                    XCAL_PROP_GAIN: 0B           , $
                    MTM_STOW: 0B                 , $
                    FAST_SLOW_CMD: 0B            , $
                    LONG_SHORT_CMD: 0B           , $
                    XCAL_STOW: 0B                , $
                    XCAL_LATCH_CMD: 0B           , $
                    MTM_LATCH_CMD: 0B            , $
                    STATMON_SPARES: bytarr( 3 )  , $
                    EXTRA_BIT_SIDE: 0B               }
 
       idx_sync = { IDX_SYNC          , $
                    SIDE: bytarr( 2 )     }
 
       ipdu_relay = { IPDU_RELAY            , $
                      POWER_RELAY4_A1: 0B   , $
                      POWER_RELAY4_B1: 0B   , $
                      POWER_RELAY4_A2: 0B   , $
                      POWER_RELAY4_B2: 0B   , $
                      POWER_RELAY2_A3: 0B   , $
                      POWER_RELAY2_B3: 0B   , $
                      OTHER_STATUS_A4: 0B   , $
                      OTHER_STATUS_B4: 0B   , $
                      IPDU_SPARES: bytarr( 4 )  }
 
       ipdu_stat = { IPDU_STAT             , $
                     MTM_POWER_ENAB: 0B    , $
                     MAIN_ELEC_CNV_H: 0B   , $
                     MAIN_ELEC_CNV_L: 0B   , $
                     HOT_SPOT_HEATER: 0B   , $
                     MTM_DRIVE_MOTOR: 0B   , $
                     XCAL_TEMP_CTRL: 0B    , $
                     SKYH_TEMP_CTRL: 0B    , $
                     REFH_TEMP_CTRL: 0B    , $
                     ICAL_TEMP_CTRL: 0B    , $
                     MTM_LATCH_MOTOR: 0B   , $
                     XCAL_LTCH_MOTOR: 0B   , $
                     XCAL_DRVE_MOTOR: 0B   , $
                     XCAL_POWER_ENAB: 0B   , $
                     LMAC_DC_CONVERT: 0B   , $
                     LMAC_AC_CONVERT: 0B   , $
                     SPARE1: 0B            , $
                     SPARE2: 0B            , $
                     SPARE3: 0B            , $
                     SPARE4: 0B                }
 
       misc_stat = { MISC_STAT                 , $
                     DWELL_FLAG_A: 0B          , $
                     DWELL_ADDRESS_A: 0B       , $
                     DWELL_FLAG_B: 0B          , $
                     DWELL_ADDRESS_B: 0B       , $
                     LVDT_STAT: bytarr( 2 )    , $
                     STATUS_RDOUT_RH: 0B       , $
                     STATUS_RDOUT_RL: 0B       , $
                     STATUS_RDOUT_LH: 0B       , $
                     STATUS_RDOUT_LL: 0B       , $
                     STAT_SPARES1: bytarr( 4 ) , $
                     INT_CONVERTER_A: 0B       , $
                     A_D_CONVERTER_A: 0B       , $
                     BIAS_PRE_CONV_A: 0B       , $
                     MTM_AND_XCAL_A: 0B        , $
                     INT_CONVERTER_B: 0B       , $
                     A_D_CONVERTER_B: 0B       , $
                     BIAS_PRE_CONV_B: 0B       , $
                     MTM_AND_XCAL_B: 0B        , $
                     STAT_SPARES2: bytarr( 8 ) , $
                     NTCH_FLTR_A: bytarr( 5 )  , $
                     NTCH_FLTR_B: bytarr( 5 )      }
 
       idx_analog = { IDX_ANALOG                   , $
                      IPDU_TEMP_A: 0B              , $
                      IPDU_TEMP_B: 0B              , $
                      DRIVEBOX_TEMP_A: 0B          , $
                      CHANNEL_PRE_AMP: 0B          , $
                      CHANNEL_TEMP_RH: 0B          , $
                      CHANNEL_TEMP_RL: 0B          , $
                      CHANNEL_TEMP_LH: 0B          , $
                      CHANNEL_TEMP_LL: 0B          , $
                      STAT_MON_TEMP_A: 0B          , $
                      STAT_MON_TEMP_B: 0B          , $
                      HOT_SPOT_CURR_A: 0B          , $
                      HOT_SPOT_CURR_B: 0B          , $
                      OPTICAL_PRE_AMP: 0B          , $
                      DRIVEBOX_TEMP_B: 0B          , $
                      ANALOG_SPARES_1: bytarr( 4 ) , $
                      IPDU_VOLTAGE_1: 0B           , $
                      IPDU_VOLTAGE_2: 0B           , $
                      IPDU_VOLTAGE_3: 0B           , $
                      IPDU_VOLTAGE_4: 0B           , $
                      IPDU_VOLTAGE_5: 0B           , $
                      IPDU_VOLTAGE_6: 0B           , $
                      IPDU_VOLTAGE_7: 0B           , $
                      IPDU_VOLTAGE_8: 0B           , $
                      IPDU_VOLTAGE_9: 0B           , $
                      IPDU_VOLTAGE_10: 0B          , $
                      IPDU_VOLTAGE_11: 0B          , $
                      IPDU_VOLTAGE_12: 0B          , $
                      IPDU_VOLTAGE_13: 0B          , $
                      IPDU_VOLTAGE_14: 0B          , $
                      IPDU_VOLTAGE_15: 0B          , $
                      IPDU_VOLTAGE_16: 0B          , $
                      IPDU_VOLTAGE_17: 0B          , $
                      IPDU_VOLTAGE_18: 0B          , $
                      IPDU_VOLTAGE_19: 0B          , $
                      IPDU_VOLTAGE_20: 0B          , $
                      IPDU_CURRENT_1: 0B           , $
                      IPDU_CURRENT_2: 0B           , $
                      IPDU_CURRENT_3: 0B           , $
                      IPDU_CURRENT_4: 0B           , $
                      IPDU_CURRENT_5: 0B           , $
                      IPDU_CURRENT_6: 0B           , $
                      IPDU_CURRENT_7: 0B           , $
                      IPDU_CURRENT_8: 0B           , $
                      IPDU_CURRENT_9: 0B           , $
                      IPDU_CURRENT_10: 0B          , $
                      IPDU_CURRENT_11: 0B          , $
                      IPDU_CURRENT_12: 0B          , $
                      ANALOG_SPARES_2: bytarr( 4 ) , $
                      MTM_CAL_MOTOR_A: 0B          , $
                      MTM_CAL_MOTOR_B: 0B              }
 
       idx_tail = { IDX_TAIL                  , $
                    CNT_TEMP_LO_TOL: 0B       , $
                    CNT_TEMP_HI_TOL: 0B       , $
                    DIH_TEMP_LO_TOL: 0B       , $
                    DIH_TEMP_HI_TOL: 0B       , $
                    BOL_TEMP_LO_TOL: 0B       , $
                    BOL_TEMP_HI_TOL: 0B       , $
                    MIR_TEMP_LO_TOL: 0B       , $
                    MIR_TEMP_HI_TOL: 0B       , $
                    TEMP_CNTRL_TOL: 0B        , $
                    OTHER_TEMP_TOL: 0B        , $
                    BIAS_RDOUT_TOL: 0B        , $
                    DIG_VOL_TOL: 0B           , $
                    ANLG_VOL_TOL: 0B          , $
                    PRE_REG_VOL_TOL: 0B       , $
                    INT_PWR_28V_TOL: 0B       , $
                    INT_PWR_15V_TOL: 0B       , $
                    INT_PWR_5V_TOL: 0B        , $
                    PRE_REG_CUR_TOL: 0B       , $
                    ANLG_CURR_TOL: 0B         , $
                    DIG_CURR_TOL: 0B          , $
                    CONST_CURR_TOL: 0B        , $
                    INT_CONV_TOL: 0B          , $
                    IDXTAIL_SPARES: bytarr( 12 )  }
 
       fdq_idx = { FDQ_IDX                              , $
                   HEADER:header                        , $
                   XCAL:xcal                            , $
                   ICHAN: replicate( ichan, 4 )         , $
                   GRTS: replicate( grts, 4 )           , $
                   TEMP_CTRL:temp_ctrl                  , $
                   STAT_MON: replicate( stat_mon, 2 )   , $
                   IDX_SYNC:idx_sync                    , $
                   IPDU_RELAY:ipdu_relay                , $
                   IPDU_STAT: replicate( ipdu_stat, 2 ) , $
                   MISC_STAT:misc_stat                  , $
                   IDX_ANALOG:idx_analog                , $
                   IDX_TAIL:idx_tail                        }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FDQ_IDX}, Nstruct )
end
