function FEX_AV_CALRS_st, Nstruct
 
  common FEX_AV_CALRS_st, defined
 
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
 
       fex_av_calrs = { FEX_AV_CALRS                      , $
                        CT_HEAD:ct_head                   , $
                        DATA_STOP: bytarr( 14 )           , $
                        DATA_STOP_TIME: Lonarr( 2 )       , $
                        PREV_DATA_START: bytarr( 14 )     , $
                        AVE_PERIOD: 0.0                   , $
                        NUM_GOOD_RECORD: 0L               , $
                        NUM_BAD_RECORD: 0L                , $
                        CALRES_AVE_A_LO: fltarr( 4 )      , $
                        CALRES_AVE_A_HI: fltarr( 4 )      , $
                        CALRES_AVE_B_LO: fltarr( 4 )      , $
                        CALRES_AVE_B_HI: fltarr( 4 )      , $
                        CALRES_DEV_A_LO: fltarr( 4 )  , $
                        CALRES_DEV_A_HI: fltarr( 4 )  , $
                        CALRES_DEV_B_LO: fltarr( 4 )  , $
                        CALRES_DEV_B_HI: fltarr( 4 )  , $
                        CALRES_BAD_A_LO: intarr( 4 ) , $
                        CALRES_BAD_A_HI: intarr( 4 ) , $
                        CALRES_BAD_B_LO: intarr( 4 ) , $
                        CALRES_BAD_B_HI: intarr( 4 ) , $
                        SPARES: bytarr( 48 )                  }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_AV_CALRS}, Nstruct )
end
