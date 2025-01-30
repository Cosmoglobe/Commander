function FEX_ENGLIM_st, Nstruct
 
  common FEX_ENGLIM_st, defined
 
  if N_elements( defined ) NE 1 then begin
 
       ct_head = { CT_HEAD                     , $
                   GMT: bytarr( 14 )           , $
                   TIME: Lonarr( 2 )           , $
                   SPACE_TIME: bytarr( 6 )     , $
                   MJR_FRM_NO: 0L              , $
                   ORBIT: 0L                   , $
                   HSKP1_TLM_FMT: 0B           , $
                   HSKP2_TLM_FMT: 0B           , $
                   INGEST_SPARES: bytarr( 18 ) , $
                   DATASET_ID: 0               , $
                   INSTR_SPARES: bytarr( 6 )       }
 
       sci_time = { SCI_TIME              , $
                    BIN_TIME: Lonarr( 2 )     }
 
       en_head = { EN_HEAD                            , $
                   LIMITS: bytarr( 23 )               , $
                   CONVERT: bytarr( 23 )              , $
                   A_SPARE: 0L                        , $
                   KEY_ID: 0                          , $
                   NUMBER_OF_RECS: 0L                 , $
                   ENGHEAD_SPARES1: bytarr( 8 )       , $
                   SCI_TIME: replicate( sci_time, 4 ) , $
                   FIRST_ENG_TIME: Lonarr( 2 )        , $
                   LAST_ENG_TIME: Lonarr( 2 )         , $
                   DQ_SUM_FLAG: bytarr( 4 )           , $
                   IFG_NO: Lonarr( 4 )                , $
                   ATT_SUM_FLAG: bytarr( 4 )          , $
                   ENGHEAD_SPARES2: bytarr( 34 )          }
 
       en_stat = { EN_STAT                       , $
                   STAT_WORD_1: 0                , $
                   INT_REF_TEMP_A: 0             , $
                   REF_HRN_TEMP_A: 0             , $
                   STAT_WORD_4: 0                , $
                   STAT_WORD_5: 0                , $
                   SKY_HRN_TEMP_A: 0             , $
                   EXT_CAL_TEMP_A: 0             , $
                   STAT_WORD_8: 0                , $
                   STAT_WORD_9: 0                , $
                   INT_REF_TEMP_B: 0             , $
                   REF_HRN_TEMP_B: 0             , $
                   STAT_WORD_12: 0               , $
                   STAT_WORD_13: 0               , $
                   SKY_HRN_TEMP_B: 0             , $
                   EXT_CAL_TEMP_B: 0             , $
                   STAT_WORD_16: 0               , $
                   GRT_ADDR: bytarr( 2 )         , $
                   MICRO_STAT_BUS: bytarr( 4 )   , $
                   BOL_CMD_BIAS: bytarr( 4 )     , $
                   DWELL_STAT: bytarr( 2 )       , $
                   LVDT_STAT: bytarr( 2 )        , $
                   HOT_SPOT_CMD: bytarr( 2 )     , $
                   ENGSTAT_SPARES1: bytarr( 10 ) , $
                   POWER_A_STATUS: bytarr( 2 )   , $
                   POWER_B_STATUS: bytarr( 2 )   , $
                   ENGSTAT_SPARES2: bytarr( 6 )      }
 
       en_analog = { EN_ANALOG                  , $
                     A_LO_GRT: fltarr( 16 )     , $
                     A_HI_GRT: fltarr( 16 )     , $
                     B_LO_GRT: fltarr( 16 )     , $
                     B_HI_GRT: fltarr( 16 )     , $
                     TEMP_CTRL: fltarr( 8 )     , $
                     IPDU_TEMP: fltarr( 2 )     , $
                     CNA_TEMP: fltarr( 4 )      , $
                     DBX_TEMP: fltarr( 2 )      , $
                     STAT_MON_TEMP: fltarr( 2 ) , $
                     PAMP_CHAN: 0.0             , $
                     PAMP_OP: 0.0               , $
                     HOT_SPOT: fltarr( 2 )      , $
                     MTM_CAL_MTR: fltarr( 2 )   , $
                     MTM_POS: fltarr( 2 )       , $
                     BOL_VOLT: fltarr( 4 )      , $
                     IPDU_VOLT: fltarr( 20 )    , $
                     IPDU_AMP: fltarr( 12 )         }
 
       en_xcal = { EN_XCAL                   , $
                   POS: intarr( 2 )          , $
                   XCAL_SPARES: bytarr( 50 )     }
 
       chan = { CHAN                       , $
                UP_SCI_MODE: 0B            , $
                FAKEIT: 0B                 , $
                UP_ADDS_GROUP: 0B          , $
                UP_SWPS_IFG: 0B            , $
                XMIT_MTM_SPEED: 0B         , $
                XMIT_MTM_LEN: 0B           , $
                SCI_GAIN: 0                , $
                DITHER: 0B                 , $
                SETUP_SPARES: bytarr( 15 )     }
 
       en_tempdiff = { EN_TEMPDIFF            , $
                       XCAL: 0                , $
                       ICAL: 0                , $
                       SKYHORN: 0             , $
                       REFHORN: 0             , $
                       DIHEDRAL: 0            , $
                       COLLIMATOR_MIRR: 0     , $
                       BOL_ASSEM: intarr( 4 )     }
 
       en_tail = { EN_TAIL                       , $
                   ENGTAIL_SPARES: bytarr( 8 )   , $
                   HSKP_FLAG: 0B                 , $
                   LMAC_ANALO_TEMP: 0.0          , $
                   LMAC_DIGIT_TEMP: 0.0          , $
                   TLM_Q_MAJ_FRM: bytarr( 2 )    , $
                   ENG_SPARES: bytarr( 9 )           }
 
       lim = { LIM                                      , $
               EN_HEAD:en_head                          , $
               EN_STAT:en_stat                          , $
               EN_ANALOG:en_analog                      , $
               EN_XCAL:en_xcal                          , $
               CHAN: replicate( chan, 4 )               , $
               EN_TEMPDIFF: replicate( en_tempdiff, 2 ) , $
               EN_TAIL:en_tail                              }
 
       fex_englim = { FEX_ENGLIM                   , $
                      CT_HEAD:ct_head              , $
                      LIM: replicate( lim, 4 )     , $
                      ENGLIM_SPARES: bytarr( 192 )     }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_ENGLIM}, Nstruct )
end
