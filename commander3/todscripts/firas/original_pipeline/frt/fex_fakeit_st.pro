function FEX_FAKEIT_st, Nstruct
 
  common FEX_FAKEIT_st, defined
 
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
 
       fex_fakeit = { FEX_FAKEIT                 , $
                      CT_HEAD:ct_head            , $
                      START_GMT: bytarr( 14 )    , $
                      START_TIME: Lonarr( 2 )    , $
                      STOP_GMT: bytarr( 14 )     , $
                      STOP_TIME: Lonarr( 2 )     , $
                      PREV_GMT: bytarr( 14 )     , $
                      FAKEIT: intarr( 2 )        , $
                      FAKEIT_CHANGE: bytarr( 2 ) , $
                      TLM_BAD_QUALITY: 0B        , $
                      FAKEIT_DATA_GAP: 0B        , $
                      END_OF_DATA: 0B            , $
                      SPARES: bytarr( 29 )           }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_FAKEIT}, Nstruct )
end
