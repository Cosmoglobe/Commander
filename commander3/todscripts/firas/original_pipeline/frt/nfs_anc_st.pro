function NFS_ANC_st, Nstruct
 
  common NFS_ANC_st, defined
 
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
 
       frame = { FRAME                                  , $
                 MIN_FRM_STATBIT: bytarr( 128 )         , $
                 SPARES: bytarr( 32 )                       }
 
       nfs_anc = { NFS_ANC                      , $
                   CT_HEAD:ct_head              , $
                   FRAME: replicate( frame, 2 ) , $
                   GMT_MJF2: Lonarr( 2 )            }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {NFS_ANC}, Nstruct )
end
