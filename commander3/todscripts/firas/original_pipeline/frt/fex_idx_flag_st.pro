function FEX_IDX_FLAG_st, Nstruct
 
  common FEX_IDX_FLAG_st, defined
 
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
 
       idx_flags = { IDX_FLAGS            , $
                     FLAGS: bytarr( 284 )     }
 
       fex_idx_flag = { FEX_IDX_FLAG              , $
                        CT_HEAD:ct_head           , $
                        IDX_FLAGS:idx_flags       , $
                        IDX_SPARES: bytarr( 164 )     }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_IDX_FLAG}, Nstruct )
end
