function NFS_HKP_st, Nstruct
 
  common NFS_HKP_st, defined
 
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
 
       hskp_head = { HSKP_HEAD                     , $
                     STAT_MONIT_CMD: intarr( 8 )   , $
                     IPDU_STAT: bytarr( 8 )        , $
                     DWELL_STAT: bytarr( 2 )       , $
                     LVDT_STAT: bytarr( 2 )        , $
                     U_PROC_STAT: bytarr( 4 )      , $
                     BOL_CMD_BIAS: bytarr( 4 )         }
 
       side_amp = { SIDE_AMP                , $
                    EX_CAL: 0               , $
                    SKY_HORN: 0             , $
                    REF_HORN: 0             , $
                    IREF_SOURCE: 0          , $
                    DIHEDRAL: 0             , $
                    BOL_ASSEM: intarr( 4 )  , $
                    MIRROR: 0               , $
                    CAL_RESIST: intarr( 4 ) , $
                    EX_CAL_SEGMENT: 0       , $
                    COLIMATOR: 0                }
 
       temps = { TEMPS                              , $
                 SIDE_AMP: replicate( side_amp, 4 ) , $
                 IPDU: bytarr( 2 )                  , $
                 DRIVE_BOX_A: 0B                    , $
                 CHAN_PRE_AMP: 0B                   , $
                 CHAN_TEMP: bytarr( 4 )             , $
                 STAT_MON: bytarr( 2 )              , $
                 HOT_SPOT_CURR: bytarr( 2 )         , $
                 OPTICAL_PREAMP: 0B                 , $
                 DRIVE_BOX_B: 0B                        }
 
       v_and_i = { V_AND_I                       , $
                   DIG_CONV_N_15: bytarr( 2 )    , $
                   DIG_CONV_P_15: bytarr( 2 )    , $
                   DIG_CONV_P_5: bytarr( 2 )     , $
                   ANA_CONV_P_15: bytarr( 2 )    , $
                   ANA_CONV_N_15: bytarr( 2 )    , $
                   BIAS_PRE_REG: bytarr( 2 )     , $
                   INT_PS_P_28: bytarr( 2 )      , $
                   INT_PS_P_15: bytarr( 2 )      , $
                   INT_PS_N_15: bytarr( 2 )      , $
                   INT_PS_P_5: bytarr( 2 )       , $
                   CUR_BIAS_PREREG: bytarr( 2 )  , $
                   CUR_ANA_CONV: bytarr( 2 )     , $
                   CUR_DIG_CONV: bytarr( 2 )     , $
                   CON_CURRENT: bytarr( 4 )      , $
                   CON_INT_CONV: bytarr( 2 )     , $
                   MTM_CAL_MOTOR: bytarr( 2 )    , $
                   BOL_BIAS_VOLT: bytarr( 4 )        }
 
       frame = { FRAME               , $
                 HSKP_HEAD:hskp_head , $
                 TEMPS:temps         , $
                 V_AND_I:v_and_i         }
 
       hskp_tail = { HSKP_TAIL              , $
                     GMT_MJF2: bytarr( 14 ) , $
                     SPARES: 0                  }
 
       mj_frm = { MJ_FRM                  , $
                  NTCH_FLT_A: bytarr( 5 ) , $
                  NTCH_FLT_B: bytarr( 5 ) , $
                  TLM_Q_MAJ_FRM: 0B       , $
                  POWER_A_STATUS: 0B      , $
                  POWER_B_STATUS: 0B      , $
                  LMAC_ANALO_TEMP: 0B     , $
                  LMAC_DIGIT_TEMP: 0B     , $
                  SPARES: bytarr( 17 )        }
 
       nfs_hkp = { NFS_HKP                        , $
                   CT_HEAD:ct_head                , $
                   FRAME: replicate( frame, 2 )   , $
                   HSKP_TAIL:hskp_tail            , $
                   MJ_FRM: replicate( mj_frm, 2 )     }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {NFS_HKP}, Nstruct )
end
