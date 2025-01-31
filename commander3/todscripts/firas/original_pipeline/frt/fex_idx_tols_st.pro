function FEX_IDX_TOLS_st, Nstruct
 
  common FEX_IDX_TOLS_st, defined
 
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
 
       idx_tols = { IDX_TOLS            , $
                    TOLS: intarr( 256 )     }
 
       fex_idx_tols = { FEX_IDX_TOLS              , $
                        CT_HEAD:ct_head           , $
                        IDX_TOLS:idx_tols         , $
                        TOL_SPARES: bytarr( 448 )     }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_IDX_TOLS}, Nstruct )
end
