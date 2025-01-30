function FEX_LIMFLAGS_st, Nstruct
 
  common FEX_LIMFLAGS_st, defined
 
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
 
       lim_flags = { LIM_FLAGS                     , $
                     FLG_BADSCI: 0B                , $
                     FLG_BADHKP: 0B                , $
                     FLG_CAL: 0B                   , $
                     FLG_MICRO: bytarr( 2 )        , $
                     FLG_ATT_ALT_RED: 0B           , $
                     FLG_ATT_ALT_YEL: 0B           , $
                     FLG_GLITCH_CT: 0B             , $
                     FLG_SATURATE: 0B              , $
                     FLG_CKSM_ERR_ST: 0B           , $
                     FLG_GRTRED_ST: bytarr( 16 )   , $
                     FLG_GRTYEL_ST: bytarr( 16 )   , $
                     FLG_TC_RED: 0B                , $
                     FLG_TC_YEL: 0B                , $
                     FLG_IPDU_TMP: bytarr( 4 )     , $
                     FLG_CHAN_TMP_ST: bytarr( 8 )  , $
                     FLG_DBOX_TMP: bytarr( 2 )     , $
                     FLG_STMON_TMP: bytarr( 4 )    , $
                     FLG_CHPA_TMP: bytarr( 2 )     , $
                     FLG_OPPA_TMP: bytarr( 2 )     , $
                     FLG_HS_HEAT_A: bytarr( 2 )    , $
                     FLG_HS_HEAT_B: bytarr( 2 )    , $
                     FLG_MTMCAL_MTR: bytarr( 4 )   , $
                     FLG_BOL_VOL_ST: bytarr( 8 )   , $
                     FLG_IPDU_VOL_ST: bytarr( 40 ) , $
                     FLG_IPDU_CUR_ST: bytarr( 24 ) , $
                     FLG_STCHG_MJ_ST: bytarr( 4 )  , $
                     FLG_DIRBE_ST: bytarr( 8 )     , $
                     FLG_DMR_ST: bytarr( 6 )       , $
                     FLG_SPACECR_ST: bytarr( 8 )   , $
                     FLG_ACLMAC_TMP: bytarr( 2 )   , $
                     FLG_DCLMAC_TMP: bytarr( 2 )   , $
                     FLG_SPARES: bytarr( 74 )      , $
                     FLG_NO_ATTITUDE: bytarr( 2 )  , $
                     FLG_ATT_SUM: bytarr( 2 )      , $
                     FLG_SUMMARY: bytarr( 2 )          }
 
       fex_limflags = { FEX_LIMFLAGS                   , $
                        CT_HEAD:ct_head                , $
                        LIM_FLAGS:lim_flags            , $
                        LIMFLAGS_SPARES: bytarr( 192 )     }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_LIMFLAGS}, Nstruct )
end
